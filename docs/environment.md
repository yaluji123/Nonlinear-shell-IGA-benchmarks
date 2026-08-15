# Runtime environment

## Tested target

The bundled plug-in files target:

- Windows x64;
- AutoCAD 2025 x64;
- ObjectARX 2025;
- Release configuration.

ObjectARX binaries are tied to the AutoCAD release. A different major AutoCAD
version normally requires rebuilding against the corresponding ObjectARX SDK.

## Loading the plug-in

1. Keep all files in `plugin/AutoCAD_2025_Release_x64/` together.
2. Add that directory to the AutoCAD trusted paths or otherwise configure
   AutoCAD security to permit loading it.
3. Load `IGAforCAD.arx` using `APPLOAD`.
4. Keep `JYH_IGAEntity.dbx`, `Clipper2Arx.crx`, and all runtime DLLs visible in
   the same directory or on `PATH`.

The Intel MKL/OpenMP, BLAS/LAPACK/ARPACK, and MinGW runtime libraries needed by
this build are included. Installing a full Intel oneAPI development environment
is not normally required merely to run the packaged plug-in.

## Build provenance

- source checkout at packaging: `54cdc577a96d04f52341052f03b9a24e6e4fa659`;
- `IGAforCAD.arx` build timestamp: 2026-08-15 09:10:45 local time;
- `IGAforCAD.arx` SHA-256:
  `1DD4ACEAC231CE046DC8689D9A6F4EED2D6A043D377451B6DCC50F100B636166`;
- `JYH_IGAEntity.dbx` SHA-256:
  `13B7FA3C09329D2D448D144DAB411B8FD8FEE49D2F8C85232BFDF3072753E31A`.

Debug symbols are excluded. The packaged implementation depends on commercial
CAD APIs and is not a standalone command-line solver.
