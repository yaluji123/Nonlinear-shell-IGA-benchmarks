# Nonlinear shell IGA benchmarks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21960376.svg)](https://doi.org/10.5281/zenodo.21960376)

This repository is a compact reproducibility package for CAD-native nonlinear
isogeometric analysis of trimmed and multipatch Kirchhoff–Love shells. It
contains the input models, selected response and field data, focused numerical
method checks, and a packaged AutoCAD/ObjectARX runtime.

The package contains seven verification configurations but five distinct
structural geometries. The two focused method checks either use a manufactured
numerical domain or reuse an existing structural model.

## Contents

- five analysis-ready AutoCAD DWG models and matching Abaqus input decks;
- the auxiliary Abaqus mode file required by the P2 initial geometry;
- IGA and Abaqus reference-point response curves for all five structures;
- compact peak-state field sources for the two stiffened-panel models;
- measured load–displacement data and the tabulated plastic hardening input for
  the perforated stiffened-panel test;
- numerical evidence and benchmark sources for the rank-revealing retained-domain
  CAS projection and the contraction-first analytical CAS–Nitsche Hessian;
- a 64-bit AutoCAD 2025/ObjectARX Release runtime with its dependent libraries.

Large Abaqus ODB files, complete nonlinear restart packages, debug symbols,
temporary solver records, and presentation graphics are intentionally excluded.
They are not required to inspect the supplied numerical evidence or to rerun the
input models.

## Repository layout

```text
data/
  benchmark_models/
    S1/ ... P2/
  experiment/
    perforated_stiffened_panel/
  method_checks/
    rank_revealing_cas/
    contraction_first_hessian/
plugin/
  AutoCAD_2025_Release_x64/
docs/
  environment.md
  model_matrix.md
  reproducibility.md
DATA_SCHEMA.md
CITATION.cff
LICENSE
```

## Structural models

| ID | Description | Distinctive feature |
|---|---|---|
| S1 | Single-patch plate | Prescribed initial geometric deformation |
| S2 | Two-patch version of S1 | Direct nonmatching Nitsche interface |
| S3 | Trimmed trapezoidal plate | Severe trimming and post-buckling response |
| P1 | Perforated stiffened panel | Explicit endplates and triangular knee braces |
| P2 | Perforated stiffened panel | Fixed and reference-point-coupled end traces |

For every model, `iga/model.dwg` is the CAD-native input and
`abaqus/model.inp` is the independent finite-element input. P2 additionally
requires `abaqus/Job-44.fil` in the Abaqus job directory because its input deck
references the first linear buckling mode stored in that file.

P1 includes explicit stiff end structures. P2 transfers their boundary action
to shell traces. One trace is fixed and the opposite trace is coupled to a
reference point. Only axial translation is released at the loaded reference
point. Compatible Kirchhoff–Love boundary-slope constraints represent the
prescribed rotations without introducing an artificial drilling degree of
freedom.

## Focused method checks

`rank_revealing_cas/` is a manufactured retained-domain test and therefore has
no separate DWG or Abaqus input. It varies the retained fraction of a fixed
quadratic candidate space and records the Gram-system conditioning, revealed
rank, finite projection status, and constant-strain reproduction error.

`contraction_first_hessian/` verifies the analytical Hessian kernel and a full
structural solve. The structural calculation reuses S2. Its two folders differ
only in whether the required contracted Hessian is assembled directly or a
pointwise full Hessian is materialized before contraction.

## Runtime

The supplied ARX/CRX/DBX files are plug-ins, not standalone executables. They
must be loaded inside a compatible 64-bit AutoCAD installation. Keep all files
in `plugin/AutoCAD_2025_Release_x64/` together. See
[`docs/environment.md`](docs/environment.md) for details.

## Integrity

`SHA256SUMS-v1.0.0.txt` contains a SHA-256 digest for every distributed file
except the checksum file itself. File meanings and column semantics are listed
in [`DATA_SCHEMA.md`](DATA_SCHEMA.md).

## Citation and use

If you use this package, cite version 1.0.0 using the version-specific DOI
[`10.5281/zenodo.21960377`](https://doi.org/10.5281/zenodo.21960377). The
concept DOI [`10.5281/zenodo.21960376`](https://doi.org/10.5281/zenodo.21960376)
resolves to the latest archived version. Complete author and repository
metadata are provided in `CITATION.cff`. The bundled plug-in is distributed
under the restricted terms in `LICENSE` because it depends on commercial CAD
APIs and third-party numerical runtimes.
