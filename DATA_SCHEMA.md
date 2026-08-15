# Data schema

## Structural model folders

Each folder under `data/benchmark_models/<ID>/` contains:

- `iga/model.dwg`: analysis-ready AutoCAD/ObjectARX model;
- `abaqus/model.inp`: Abaqus/Standard input deck;
- `response/iga_reference_point_history.csv`: accepted IGA reference-point
  states;
- `response/abaqus_reference_point_curve.csv`: Abaqus reference-point response
  extracted frame by frame;
- `response/iga_result_manifest.json`: solver and result-package metadata for
  the supplied IGA curve.

P2 also contains `abaqus/Job-44.fil`, which is referenced by its `*IMPERFECTION`
card. The filename is retained because Abaqus resolves the reference by name.

### IGA reference-point history

The CSV contains one row per accepted state. Important columns are:

- `frame_id`, `accepted_step`, `attempt`, `iterations`: state identifiers;
- `load_factor`, `pseudo_time`: incremental solution coordinates;
- `reference_point_id` and `reference_x/y/z`: reference-point identity and
  position in millimetres;
- `u1/u2/u3`: translations in millimetres;
- `ur1/ur2/ur3`: rotations in radians;
- `rf1/rf2/rf3`: reactions in newtons;
- `rm1/rm2/rm3`: reaction moments in newton millimetres.

### Abaqus reference-point curve

The CSV contains one row per output frame. `displacement_mm` and `reaction_n`
are the positive magnitudes used for response-curve comparison. The signed
components and all six reference-point reactions remain available in the
remaining columns.

## Peak-state field sources

P1 and P2 contain compact source data under `peak_state/`.

### IGA

- `mesh.bin`: nonlinear result mesh;
- `frame_XXXXXXXX.bin`: the accepted state identified as the peak-response
  frame in `source_result_manifest.json`;
- `source_result_manifest.json`: frame numbering, field semantics, units, and
  peak bracketing metadata.

Only the required peak frame is distributed; this is not a complete restartable
result package.

### Abaqus

- `nodes.csv`: original coordinates, displacement components and deformed
  coordinates in millimetres;
- `elements.csv`: element connectivity;
- `nodal_mises_display.csv`: displayed nodal Mises stress in MPa, obtained by
  taking the maximum across section points at each element node and then
  averaging contributions from adjacent elements.

The displayed nodal stress is intended for contour comparison. It is not a
replacement for integration-point stress when checking constitutive response.

## Experimental data

`data/experiment/perforated_stiffened_panel/load_displacement.xlsx` contains
the measured axial load and end-shortening records. The workbook is preserved
in its supplied form.

`material_hardening.csv` contains the true yield stress in MPa and equivalent
plastic strain used by the P2 numerical model. The same values are present in
the Abaqus input deck.

## Rank-revealing retained-domain CAS check

`rank_revelation.csv` records:

- retained fraction;
- full Gram-system status and condition number;
- retained numerical rank selected by the production projection;
- finite projection status;
- weighted constant-strain reproduction error.

`rank_revelation.json` contains the underlying benchmark record.

## Contraction-first Hessian check

`hessian_scaling.csv` compares full pointwise Hessian materialization with
direct contraction as the equivalent workload is increased. It reports memory,
median kernel time, maximum component difference, and equivalent-physics
difference.

`cad_hessian_equivalence.csv` reports wall time, interface-kernel time, active
pointwise columns, and analytical publication batches for the two structural
executions. `cad_hessian_comparison.json` contains the detailed path and
response checks. The `structural_runs/` folders preserve the accepted histories,
reference-point curves, manifests, and runtime records.
