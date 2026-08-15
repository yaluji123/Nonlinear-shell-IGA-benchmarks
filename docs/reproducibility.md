# Reproducibility notes

## Inspection without AutoCAD

The CSV, JSON, XLSX, INP, and C++ benchmark files can be inspected without the
CAD runtime. They support independent checks of:

- load–displacement response histories;
- peak-response state selection;
- retained-domain rank detection and constant-strain reproduction;
- Hessian component equivalence, temporary storage, and kernel timing;
- measured plastic hardening and physical load–displacement data.

## CAD-native rerun

The DWG files are analysis-ready inputs for the packaged plug-in. Load the
runtime inside AutoCAD 2025 x64, open a copied DWG, and run the nonlinear solver
from the plug-in interface. Use a fresh working copy because AutoCAD may update
DWG metadata when a drawing is opened or saved.

The model definitions already store their geometry, shell properties,
constraints, reference points, interface entities, material tables, and
nonlinear solution settings. S1–S3 and P1–P2 may therefore be run through the
same user interface used to inspect them.

## Abaqus rerun

Run each `abaqus/model.inp` with Abaqus/Standard. P2 must be run in a directory
that also contains `Job-44.fil`; its input deck references that exact basename.
The Abaqus inputs are independent comparison models and are not consumed by the
AutoCAD plug-in.

## Focused method checks

The two C++ files under `data/method_checks/` call the same numerical kernels
used to form the supplied evidence. They are provided for inspection and
rebuilding in an environment that has the required compiler and numerical
dependencies. The committed CSV and JSON files are the frozen outputs of the
distributed configuration.

## Exclusions

Complete ODB databases and restartable nonlinear result directories are not
included because they are large generated artifacts. The distributed response
curves and compact peak-state sources contain the quantities needed for the
documented comparisons. Full result files can be regenerated from the supplied
inputs in licensed AutoCAD and Abaqus environments.
