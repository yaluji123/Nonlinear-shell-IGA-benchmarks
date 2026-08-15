#include "../IGAforCAD/TrimAwareSymmetricCAS.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kBasisCount = 9u;
constexpr double kRelativeRankTolerance = 1.0e-11;

const char* StatusLabel(iga::nonlinear::AssemblyStatus status)
{
	switch (status) {
	case iga::nonlinear::AssemblyStatus::Ok:
		return "Ok";
	case iga::nonlinear::AssemblyStatus::InvalidState:
		return "InvalidState";
	case iga::nonlinear::AssemblyStatus::ElementFailure:
		return "ElementFailure";
	}
	return "Unknown";
}

std::array<double, 8> GaussPoints()
{
	return { -0.9602898564975363, -0.7966664774136267,
		-0.5255324099163290, -0.1834346424956498,
		0.1834346424956498, 0.5255324099163290,
		0.7966664774136267, 0.9602898564975363 };
}

std::array<double, 8> GaussWeights()
{
	return { 0.1012285362903763, 0.2223810344533745,
		0.3137066458778873, 0.3626837833783620,
		0.3626837833783620, 0.3137066458778873,
		0.2223810344533745, 0.1012285362903763 };
}

std::array<double, kBasisCount> TensorQuadraticBasis(double x, double y)
{
	const std::array<double, 3> px = { 1.0, x, x * x };
	const std::array<double, 3> py = { 1.0, y, y * y };
	std::array<double, kBasisCount> values{};
	for (std::size_t i = 0u; i < 3u; ++i) {
		for (std::size_t j = 0u; j < 3u; ++j)
			values[3u * i + j] = px[i] * py[j];
	}
	return values;
}

std::array<double, kBasisCount> SymmetricEigenvalues(
	std::array<double, kBasisCount * kBasisCount> matrix)
{
	for (int sweep = 0; sweep < 120; ++sweep) {
		std::size_t p = 0u;
		std::size_t q = 1u;
		double maximum = 0.0;
		for (std::size_t i = 0u; i < kBasisCount; ++i) {
			for (std::size_t j = i + 1u; j < kBasisCount; ++j) {
				const double value = std::abs(matrix[i * kBasisCount + j]);
				if (value > maximum) {
					maximum = value;
					p = i;
					q = j;
				}
			}
		}
		if (maximum <= 1.0e-28)
			break;
		const double app = matrix[p * kBasisCount + p];
		const double aqq = matrix[q * kBasisCount + q];
		const double apq = matrix[p * kBasisCount + q];
		const double angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
		const double c = std::cos(angle);
		const double s = std::sin(angle);
		for (std::size_t k = 0u; k < kBasisCount; ++k) {
			if (k == p || k == q)
				continue;
			const double mkp = matrix[k * kBasisCount + p];
			const double mkq = matrix[k * kBasisCount + q];
			matrix[k * kBasisCount + p] = c * mkp - s * mkq;
			matrix[p * kBasisCount + k] = matrix[k * kBasisCount + p];
			matrix[k * kBasisCount + q] = s * mkp + c * mkq;
			matrix[q * kBasisCount + k] = matrix[k * kBasisCount + q];
		}
		matrix[p * kBasisCount + p] = c * c * app
			- 2.0 * s * c * apq + s * s * aqq;
		matrix[q * kBasisCount + q] = s * s * app
			+ 2.0 * s * c * apq + c * c * aqq;
		matrix[p * kBasisCount + q] = 0.0;
		matrix[q * kBasisCount + p] = 0.0;
	}
	std::array<double, kBasisCount> eigenvalues{};
	for (std::size_t i = 0u; i < kBasisCount; ++i)
		eigenvalues[i] = matrix[i * kBasisCount + i];
	std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<double>());
	return eigenvalues;
}

struct CaseResult
{
	double retainedFraction = 0.0;
	const char* productionStatus = "Unknown";
	std::size_t retainedRank = 0u;
	double conditionNumber = std::numeric_limits<double>::infinity();
	double constantError = std::numeric_limits<double>::infinity();
	bool finite = false;
	const char* fullRankStatus = "UNSTABLE";
};

CaseResult EvaluateCase(double retainedFraction)
{
	const auto points = GaussPoints();
	const auto gaussWeights = GaussWeights();
	std::vector<double> physicalWeights;
	std::vector<double> assumedBasisValues;
	std::vector<double> compatibleValues;
	physicalWeights.reserve(64u);
	assumedBasisValues.reserve(64u * kBasisCount);
	compatibleValues.reserve(64u);
	std::array<double, kBasisCount * kBasisCount> gram{};

	// The retained physical strip is x in [-1,-1+2*rho], y in [-1,1].
	// This is the canonical cut-cell sliver: the parent polynomial basis is
	// regular, but its restriction becomes nearly rank deficient as rho -> 0.
	for (std::size_t i = 0u; i < points.size(); ++i) {
		const double x = -1.0 + retainedFraction * (points[i] + 1.0);
		for (std::size_t j = 0u; j < points.size(); ++j) {
			const double y = points[j];
			const double weight = retainedFraction
				* gaussWeights[i] * gaussWeights[j];
			const auto basis = TensorQuadraticBasis(x, y);
			physicalWeights.push_back(weight);
			assumedBasisValues.insert(
				assumedBasisValues.end(), basis.begin(), basis.end());
			compatibleValues.push_back(1.0);
			for (std::size_t a = 0u; a < kBasisCount; ++a) {
				for (std::size_t b = 0u; b < kBasisCount; ++b)
					gram[a * kBasisCount + b] += weight * basis[a] * basis[b];
			}
		}
	}

	iga::nonlinear::SymmetricCASProjectionInput input;
	input.sampleCount = physicalWeights.size();
	input.assumedBasisCount = kBasisCount;
	input.valueCount = 1u;
	input.physicalWeights = physicalWeights.data();
	input.assumedBasisValues = assumedBasisValues.data();
	input.compatibleValues = compatibleValues.data();
	input.relativeRankTolerance = kRelativeRankTolerance;
	iga::nonlinear::SymmetricCASProjectionOutput output;
	const auto status = iga::nonlinear::ProjectTrimAwareSymmetricCAS(input, output);

	double errorSquared = 0.0;
	bool finite = status == iga::nonlinear::AssemblyStatus::Ok;
	if (finite) {
		for (std::size_t sample = 0u; sample < physicalWeights.size(); ++sample) {
			const double value = output.projectedValues[sample];
			finite = finite && std::isfinite(value);
			const double difference = value - 1.0;
			errorSquared += physicalWeights[sample] * difference * difference;
		}
	}

	const auto eigenvalues = SymmetricEigenvalues(gram);
	const double largest = std::max(0.0, eigenvalues.front());
	const double smallest = std::max(0.0, eigenvalues.back());
	const double condition = smallest > 0.0
		? largest / smallest
		: std::numeric_limits<double>::infinity();
	const bool stableFullRank = std::isfinite(condition) && condition < 1.0e14;

	CaseResult result;
	result.retainedFraction = retainedFraction;
	result.productionStatus = StatusLabel(status);
	result.retainedRank = output.retainedRank;
	result.conditionNumber = condition;
	result.constantError = std::sqrt(std::max(0.0, errorSquared));
	result.finite = finite;
	result.fullRankStatus = stableFullRank ? "STABLE" : "ILL_CONDITIONED";
	return result;
}

void WriteNumber(std::ostream& stream, double value)
{
	if (std::isfinite(value))
		stream << std::setprecision(17) << value;
	else
		stream << "1e999";
}

} // namespace

int main(int argc, char** argv)
{
	if (argc != 3 || std::string(argv[1]) != "--output") {
		std::cerr << "usage: benchmark_trim_aware_cas_rank_revelation --output FILE\n";
		return 2;
	}
	const std::array<double, 4> retainedFractions = {
		0.5, 0.01, 0.0001, 0.000001 };
	std::array<CaseResult, retainedFractions.size()> results{};
	for (std::size_t i = 0u; i < retainedFractions.size(); ++i)
		results[i] = EvaluateCase(retainedFractions[i]);
	const bool severeRequired =
		results.back().productionStatus == std::string("Ok")
		&& results.back().retainedRank < kBasisCount
		&& results.back().fullRankStatus != std::string("STABLE");

	std::ofstream output(argv[2], std::ios::binary | std::ios::trunc);
	if (!output) {
		std::cerr << "could not open output\n";
		return 3;
	}
	output << "{\n"
		<< "  \"schema\": \"iga.trim-cas-rank-revelation-benchmark.v1\",\n"
		<< "  \"production_seam\": \"ProjectTrimAwareSymmetricCAS\",\n"
		<< "  \"assumed_basis_count\": " << kBasisCount << ",\n"
		<< "  \"relative_rank_tolerance\": "
		<< std::setprecision(17) << kRelativeRankTolerance << ",\n"
		<< "  \"cases\": [\n";
	for (std::size_t i = 0u; i < results.size(); ++i) {
		const CaseResult& result = results[i];
		output << "    {\"retained_fraction\": ";
		WriteNumber(output, result.retainedFraction);
		output << ", \"production_status\": \"" << result.productionStatus
			<< "\", \"retained_rank\": " << result.retainedRank
			<< ", \"gram_condition_number\": ";
		WriteNumber(output, result.conditionNumber);
		output << ", \"weighted_constant_reproduction_error\": ";
		WriteNumber(output, result.constantError);
		output << ", \"projection_values_finite\": "
			<< (result.finite ? "true" : "false")
			<< ", \"full_rank_solve_status\": \""
			<< result.fullRankStatus << "\"}"
			<< (i + 1u == results.size() ? "\n" : ",\n");
	}
	output << "  ],\n"
		<< "  \"conclusion\": {\"severe_sliver_rank_revelation_required\": "
		<< (severeRequired ? "true" : "false") << "}\n"
		<< "}\n";
	return severeRequired ? 0 : 4;
}
