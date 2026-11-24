# Validation status

The standalone RC-TLR module has been statically inspected. MATLAB R2024b
is not available in the preparation environment, so neither the MATLAB
tests nor the full 55-case experimental suite have been executed here.
Run `test_biomol_rc_tlr()` in MATLAB before collecting benchmark data.

## Checks performed

- Matched the leading MATLAB function name to every source filename and
  scanned for legacy approximation dependencies and unresolved old method names.
- Independently counted the deduplicated default plan: 55 cases, using only
  the CP-input and full-grid RC-TLR routes (including the optional lattice).
- Loaded the nine original biomolecular data files with SciPy and checked
  their CP factor dimensions, finite values, finite weights, and nonzero
  grid scale. The lattice input has 50 original terms, 28 retained terms
  and three length-960 grids.
- Computed each original data file's SHA-256 and verified the manifest.
  Each uploaded binary blob, including all 29 binary parts of the seven
  larger files, matched the corresponding local Git object SHA.
- Independently evaluated the weighted interpolation equations in Python
  on `n256_Pn500.mat`, at p=33, RRQR tolerance 4e-8 and strength 2.
  The selected ranks were (29,26,29), maximum weighted unfolding residual
  1.10e-8, maximum strength statistic 1.488, and unweighted Chebyshev sample
  relative Frobenius compression error 2.14e-8. Cardinal rows were exactly
  identity. This checks the equations, not execution of the MATLAB driver
  or the final uniform-grid slice error.

The included MATLAB tests cover the functional and full-grid spline
algorithms, rectangular grids, cardinal interpolation, weighted CP
orthogonal reduction, an end-to-end synthetic biomolecular input, storage
and timing identities, and checkpoint resumption. The data restoration
routine also checks byte counts and SHA-256 of original/restored MAT files.

No numerical benchmark tables, MATLAB timings or measured peak-memory
claims are made for this refactored module.
