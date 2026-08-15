# Model matrix

## Five structural geometries

### S1: single-patch plate

S1 is a single Kirchhoff–Love surface with an imposed initial geometric
deformation. It is loaded by prescribed end displacement and provides the
single-patch nonlinear reference.

### S2: directly coupled two-patch plate

S2 is formed by splitting the S1 surface into two patches while leaving the
physical geometry, material, initial geometry, loading, and end conditions
unchanged. The new common line is coupled by the direct shell Nitsche operator.

### S3: trimmed trapezoidal plate

S3 introduces a severely trimmed boundary and follows its nonlinear response
under axial shortening.

### P1: perforated stiffened panel with explicit end structures

P1 includes a perforated plate, four longitudinal stiffeners, endplates, and
triangular knee braces. The end structures create stiff zones that restrain end
rotation. Its initial geometry is the supplied three-component deformation.

### P2: perforated stiffened panel with boundary idealization

P2 uses the same perforated stiffened test section without discretizing the
massive end structures. One shell trace is fixed. The other is coupled to a
reference point with only axial translation released. The model uses a scaled
linear buckling mode as its initial geometry and a measured tabulated plastic
hardening law.

## Two focused numerical configurations

The rank-revealing CAS check is a manufactured retained-domain calculation and
does not define another structural geometry. The contraction-first Hessian
check reuses S2 and compares two mathematically equivalent assembly orders.
