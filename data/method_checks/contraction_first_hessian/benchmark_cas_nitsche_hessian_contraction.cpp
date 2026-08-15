#include "../IGAforCAD/AggregatedCASNitscheHessian.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using iga::nonlinear::AggregateCASContractionWeights;
using iga::nonlinear::AssemblyStatus;
using iga::nonlinear::CASAggregatedNitscheHessianOutput;
using iga::nonlinear::CASAggregatedNitscheHessianTransactionOptions;
using iga::nonlinear::CASBatchPatchDofLayout;
using iga::nonlinear::CASBatchSampleProjectedResultantWeights;
using iga::nonlinear::CASContractionWeightAggregationOutput;
using iga::nonlinear::CASContractedBatchHessianEvaluation;
using iga::nonlinear::CASInterfaceSlotPenaltyBlockContraction;
using iga::nonlinear::EvaluateAndScatterAggregatedCASNitscheHessians;
using iga::nonlinear::FixedElasticCASBatchAnalyticDisposition;

constexpr std::size_t kBaseSampleCount = 64u;
constexpr std::size_t kSupportControlPoints = 16u;
constexpr int kDofPerControlPoint = 3;
constexpr std::size_t kDofCount =
	kSupportControlPoints * static_cast<std::size_t>(kDofPerControlPoint);
constexpr std::size_t kMatrixEntries = kDofCount * kDofCount;

volatile double g_digest = 0.0;

double AnalyticGeneralizedForceHessian(
	std::size_t originalSample,
	std::size_t resultant,
	std::size_t row,
	std::size_t column)
{
	const std::uint64_t symmetricIndex = static_cast<std::uint64_t>(
		(std::min)(row, column) * kDofCount + (std::max)(row, column));
	const std::uint64_t mixed =
		(static_cast<std::uint64_t>(originalSample + 1u) * 1315423911ULL)
		^ (static_cast<std::uint64_t>(resultant + 3u) * 2654435761ULL)
		^ (symmetricIndex * 2246822519ULL);
	const double oscillatory = static_cast<double>(mixed % 1009ULL) / 1009.0;
	return (row == column ? 0.01 * static_cast<double>(resultant + 1u) : 0.0)
		+ 1.0e-4 * (oscillatory - 0.5);
}

std::vector<CASInterfaceSlotPenaltyBlockContraction> BuildContributions(
	std::size_t multiplier)
{
	const std::size_t sampleCount = kBaseSampleCount * multiplier;
	std::vector<CASInterfaceSlotPenaltyBlockContraction> contributions;
	contributions.reserve(sampleCount);
	for (std::size_t sample = 0u; sample < sampleCount; ++sample) {
		CASInterfaceSlotPenaltyBlockContraction contribution;
		contribution.interfaceSlotId = static_cast<std::uint64_t>(sample + 1u);
		contribution.penaltyBlockOrdinal = 0u;
		contribution.scale = 1.0 / static_cast<double>(multiplier);
		contribution.gamma1 = 1.0;
		contribution.gamma2 = 0.0;
		CASBatchSampleProjectedResultantWeights term;
		term.batchId = 101u;
		term.sampleOrdinal = sample;
		term.sampleCount = sampleCount;
		const std::size_t original = sample / multiplier;
		for (std::size_t resultant = 0u; resultant < 6u; ++resultant) {
			term.localResultantWeightsAfterProjectedRotationPullback[resultant] =
				(1.0 + static_cast<double>((original + 2u * resultant) % 7u))
				/ 17.0;
		}
		contribution.side1.push_back(term);
		contributions.push_back(std::move(contribution));
	}
	return contributions;
}

std::vector<double> FullPointwiseThenContract(
	const std::vector<double>& weights,
	std::size_t multiplier)
{
	const std::size_t sampleCount = weights.size() / 6u;
	std::vector<double> fullHessians(sampleCount * 6u * kMatrixEntries, 0.0);
	for (std::size_t sample = 0u; sample < sampleCount; ++sample) {
		const std::size_t original = sample / multiplier;
		for (std::size_t resultant = 0u; resultant < 6u; ++resultant) {
			const std::size_t offset =
				(sample * 6u + resultant) * kMatrixEntries;
			for (std::size_t row = 0u; row < kDofCount; ++row) {
				for (std::size_t column = 0u; column < kDofCount; ++column) {
					fullHessians[offset + row * kDofCount + column] =
						AnalyticGeneralizedForceHessian(
							original, resultant, row, column);
				}
			}
		}
	}
	std::vector<double> contracted(kMatrixEntries, 0.0);
	for (std::size_t sample = 0u; sample < sampleCount; ++sample) {
		for (std::size_t resultant = 0u; resultant < 6u; ++resultant) {
			const double weight = weights[sample * 6u + resultant];
			const std::size_t offset =
				(sample * 6u + resultant) * kMatrixEntries;
			for (std::size_t entry = 0u; entry < kMatrixEntries; ++entry)
				contracted[entry] += weight * fullHessians[offset + entry];
		}
	}
	return contracted;
}

std::vector<double> ContractFirst(
	const std::vector<double>& weights,
	std::size_t multiplier)
{
	const std::size_t sampleCount = weights.size() / 6u;
	std::vector<double> contracted(kMatrixEntries, 0.0);
	// Contraction-first ordering keeps only one scalar accumulator live for a
	// tangent entry and publishes that entry once.  The pointwise comparator
	// above necessarily materializes and rereads S*6*D^2 values.
	for (std::size_t row = 0u; row < kDofCount; ++row) {
		for (std::size_t column = 0u; column < kDofCount; ++column) {
			double value = 0.0;
			for (std::size_t sample = 0u; sample < sampleCount; ++sample) {
				const std::size_t original = sample / multiplier;
				for (std::size_t resultant = 0u; resultant < 6u; ++resultant) {
					value += weights[sample * 6u + resultant]
						* AnalyticGeneralizedForceHessian(
							original, resultant, row, column);
				}
			}
			contracted[row * kDofCount + column] = value;
		}
	}
	return contracted;
}

double MaximumDifference(
	const std::vector<double>& first,
	const std::vector<double>& second)
{
	if (first.size() != second.size())
		return std::numeric_limits<double>::infinity();
	double maximum = 0.0;
	for (std::size_t i = 0u; i < first.size(); ++i)
		maximum = (std::max)(maximum, std::abs(first[i] - second[i]));
	return maximum;
}

double MedianMilliseconds(std::vector<double> samples)
{
	std::sort(samples.begin(), samples.end());
	return samples[samples.size() / 2u];
}

struct BenchmarkCase
{
	std::size_t multiplier = 0u;
	double maximumDifference = 0.0;
	double physicsDifference = 0.0;
	std::size_t fullBytes = 0u;
	std::size_t contractedBytes = 0u;
	double pointwiseMilliseconds = 0.0;
	double contractionMilliseconds = 0.0;
	double speedup = 0.0;
	double memoryReduction = 0.0;
};

} // namespace

int main(int argc, char** argv)
{
	std::string outputPath;
	int repeats = 5;
	for (int index = 1; index < argc; ++index) {
		const std::string argument(argv[index]);
		if (argument == "--output" && index + 1 < argc)
			outputPath = argv[++index];
		else if (argument == "--repeats" && index + 1 < argc)
			repeats = std::stoi(argv[++index]);
		else {
			std::cerr << "invalid arguments\n";
			return 2;
		}
	}
	if (outputPath.empty() || repeats < 3 || repeats > 25) {
		std::cerr << "usage: benchmark --output FILE --repeats N(3..25)\n";
		return 2;
	}

	const std::array<std::size_t, 4> multipliers = { 1u, 2u, 4u, 8u };
	std::array<BenchmarkCase, multipliers.size()> cases{};
	std::vector<double> referencePhysics;
	std::array<int, kDofCount> localToGlobal{};
	for (std::size_t dof = 0u; dof < kDofCount; ++dof)
		localToGlobal[dof] = static_cast<int>(dof);

	for (std::size_t caseIndex = 0u; caseIndex < multipliers.size(); ++caseIndex) {
		const std::size_t multiplier = multipliers[caseIndex];
		auto contributions = BuildContributions(multiplier);
		CASContractionWeightAggregationOutput aggregation;
		if (AggregateCASContractionWeights(
				contributions.data(), contributions.size(), aggregation)
			!= AssemblyStatus::Ok || aggregation.batches.size() != 1u) {
			std::cerr << "aggregation failed\n";
			return 3;
		}
		const std::vector<double>& weights =
			aggregation.batches.front().generalizedResultantWeights;

		std::vector<double> pointwiseTimes;
		std::vector<double> contractionTimes;
		std::vector<double> pointwise;
		std::vector<double> contracted;
		for (int repeat = 0; repeat < repeats + 1; ++repeat) {
			auto started = std::chrono::steady_clock::now();
			pointwise = FullPointwiseThenContract(weights, multiplier);
			auto stopped = std::chrono::steady_clock::now();
			const double pointwiseMs = std::chrono::duration<double, std::milli>(
				stopped - started).count();

			started = std::chrono::steady_clock::now();
			contracted = ContractFirst(weights, multiplier);
			stopped = std::chrono::steady_clock::now();
			const double contractedMs = std::chrono::duration<double, std::milli>(
				stopped - started).count();
			g_digest += pointwise.front() + contracted.back();
			if (repeat != 0) {
				pointwiseTimes.push_back(pointwiseMs);
				contractionTimes.push_back(contractedMs);
			}
		}

		CASAggregatedNitscheHessianOutput publicOutput;
		auto resolve = [&](std::uint64_t batchId, CASBatchPatchDofLayout& layout) {
			if (batchId != 101u)
				return AssemblyStatus::InvalidState;
			layout.controlPointCount = kSupportControlPoints;
			layout.dofPerControlPoint = kDofPerControlPoint;
			layout.patchLocalToGlobalDofs = localToGlobal.data();
			layout.patchLocalToGlobalDofCount = localToGlobal.size();
			return AssemblyStatus::Ok;
		};
		auto evaluate = [&](std::uint64_t batchId, const double* rawWeights,
			std::size_t weightCount, CASContractedBatchHessianEvaluation& evaluation) {
			if (batchId != 101u || weightCount != weights.size())
				return AssemblyStatus::InvalidState;
			std::vector<double> copied(rawWeights, rawWeights + weightCount);
			evaluation.contractedGeneralizedForceHessian =
				ContractFirst(copied, multiplier);
			evaluation.supportControlPoints.resize(kSupportControlPoints);
			for (std::size_t cp = 0u; cp < kSupportControlPoints; ++cp)
				evaluation.supportControlPoints[cp] = static_cast<int>(cp);
			evaluation.dofPerControlPoint = kDofPerControlPoint;
			evaluation.maximumRelativeAsymmetry = 0.0;
			evaluation.analyticDisposition =
				FixedElasticCASBatchAnalyticDisposition::Published;
			return AssemblyStatus::Ok;
		};
		CASAggregatedNitscheHessianTransactionOptions options;
		options.emissionTolerance = 0.0;
		if (EvaluateAndScatterAggregatedCASNitscheHessians(
				aggregation, resolve, evaluate, options, publicOutput)
			!= AssemblyStatus::Ok
			|| publicOutput.contractedHessianValueBytes
				!= kMatrixEntries * sizeof(double)) {
			std::cerr << "public contraction-first seam failed\n";
			return 4;
		}

		if (caseIndex == 0u)
			referencePhysics = contracted;
		BenchmarkCase& result = cases[caseIndex];
		result.multiplier = multiplier;
		result.maximumDifference = MaximumDifference(pointwise, contracted);
		result.physicsDifference = MaximumDifference(referencePhysics, contracted);
		result.fullBytes = kBaseSampleCount * multiplier * 6u
			* kMatrixEntries * sizeof(double);
		result.contractedBytes = publicOutput.contractedHessianValueBytes;
		result.pointwiseMilliseconds = MedianMilliseconds(pointwiseTimes);
		result.contractionMilliseconds = MedianMilliseconds(contractionTimes);
		result.speedup = result.pointwiseMilliseconds
			/ result.contractionMilliseconds;
		result.memoryReduction = static_cast<double>(result.fullBytes)
			/ static_cast<double>(result.contractedBytes);
	}

	const bool exact = std::all_of(cases.begin(), cases.end(),
		[](const BenchmarkCase& value) {
			return value.maximumDifference <= 1.0e-11
				&& value.physicsDifference <= 1.0e-11;
		});
	const bool scaling = cases.back().speedup >= 1.25
		&& cases.back().memoryReduction >= 1000.0;
	std::ofstream output(outputPath, std::ios::binary | std::ios::trunc);
	if (!output)
		return 5;
	output << std::setprecision(17)
		<< "{\n"
		<< "  \"schema\": \"iga.cas-nitsche-contraction-benchmark.v1\",\n"
		<< "  \"production_seams\": [\"AggregateCASContractionWeights\", "
			"\"EvaluateAndScatterAggregatedCASNitscheHessians\"],\n"
		<< "  \"base_sample_count\": " << kBaseSampleCount << ",\n"
		<< "  \"support_dof_count\": " << kDofCount << ",\n"
		<< "  \"repeats\": " << repeats << ",\n"
		<< "  \"cases\": [\n";
	for (std::size_t index = 0u; index < cases.size(); ++index) {
		const BenchmarkCase& value = cases[index];
		output << "    {\"equivalent_workload_multiplier\": " << value.multiplier
			<< ", \"maximum_absolute_difference\": " << value.maximumDifference
			<< ", \"equivalent_physics_difference_from_m1\": " << value.physicsDifference
			<< ", \"pointwise_full_hessian_bytes\": " << value.fullBytes
			<< ", \"contracted_hessian_bytes\": " << value.contractedBytes
			<< ", \"memory_reduction_ratio\": " << value.memoryReduction
			<< ", \"pointwise_median_ms\": " << value.pointwiseMilliseconds
			<< ", \"contraction_first_median_ms\": " << value.contractionMilliseconds
			<< ", \"median_speedup\": " << value.speedup << "}"
			<< (index + 1u == cases.size() ? "\n" : ",\n");
	}
	output << "  ],\n"
		<< "  \"conclusion\": {\"contraction_first_exact\": "
		<< (exact ? "true" : "false")
		<< ", \"workload_scaling_supported\": "
		<< (scaling ? "true" : "false") << "}\n"
		<< "}\n";
	return exact && scaling && std::isfinite(g_digest) ? 0 : 6;
}
