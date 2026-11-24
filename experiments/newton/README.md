# Coulomb / Newton kernel experiments

This module contains the supplied MATLAB implementations for the cell-averaged Coulomb (Newton) kernel on ([-1,1]^3), with both the full canonical tensor and its long-range component. It also contains a standalone RC-TLR CP-input experiment. The original MATLAB source files are preserved without numerical or naming changes.

## Requirements

- MATLAB (R2024b is the development environment).
- The Newton canonical factor data from [ChebTuck](https://github.com/bonans/ChebTuck/tree/2aa3729fe9991377bf1a6815ab3edd73423689dd/data/Newton_CP_-1_1).
- The supplied comparison driver uses MATLAB built-in numerical routines; it does not call Chebfun or Tensor Toolbox.

## Data

Download the original `Newt_canon_<n>.mat` files from the pinned [ChebTuck dataset directory](https://github.com/bonans/ChebTuck/tree/2aa3729fe9991377bf1a6815ab3edd73423689dd/data/Newton_CP_-1_1). Each file contains the canonical factor matrix `U1r`. Keep these files under `experiments/newton/data/Newton_CP_-1_1/`, or pass the absolute path to your existing dataset directory.

Available grid sizes: 256, 512, 1024, 2048, 4096, 8192. The benchmark scales `U1r` by (1/h), where (h=2/n), and uses the cell-centred grid. The long-range part retains the first `ceil(R/2)-3` canonical columns, with (R) the full CP rank. These are cell-averaged kernel data, not point samples at the singularity.

## Quick start

From the repository root in MATLAB:

```matlab
newtonDir = fullfile(pwd,'experiments','newton');
addpath(newtonDir);
dataDir = fullfile(newtonDir,'data','Newton_CP_-1_1');

% Small comparison: both full and long-range parts, n=256, p=129.
[results, firstPass] = compare_newton_paper(dataDir,256,129,129,1);
disp(results);
disp(firstPass);
```

To run the supplied full comparison sweep:

```matlab
[results, firstPass] = compare_newton_paper(dataDir);
writetable(results,'newton_comparison.csv');
writetable(firstPass,'newton_first_passing_p.csv');
```

For the standalone RC-TLR CP-input experiment:

```matlab
opts = struct('n',256,'pValues',129,'dataDir',dataDir, ...
              'warmups',0,'repeats',1,'qrTol',5e-8);
[results,model,diagnostics] = rc_tlr_1d_experiment(opts);
disp(results);
```

Increase node counts only after checking available memory. The general builder constructs a pairwise CP-factor product with `p^2` rows, whereas the symmetric comparison driver works in reduced orthonormal coordinates. Avoid dense `kernel_tensor(n)` calls for large grids; request `'cp'` or `'slice'` instead.

## Supplied numerical settings

| Setting | Comparison driver | Standalone RC-TLR |
|---|---|---|
| Grid sizes | 256, 512, 1024, 2048, 4096 | 256 |
| Full kernel node counts | 129, 257, 513, 1025, 2049, 4097, 8193, 16385 | Not the default input |
| Long-range node counts | 129, 257, 513, 1025, 2049, 4097 | 129, 257, 513, 1025 |
| RC-TLR RRQR tolerance | 1e-9 | 5e-8 (configurable) |
| ChebTuck compression tolerance | 1e-7 | Not applicable |
| Target relative error | 1e-7 | 1e-7 |
| Strong RRQR strength | 2 | 2 (configurable) |
| Validation | Relative maximum error, complete middle slice | Same by default; full-grid Frobenius error optional |

The RRQR tolerance and the target final approximation error are distinct parameters. These defaults reflect the uploaded implementation; confirm them against the final manuscript before reproducing publication tables. In particular, the two entry points still use different RRQR defaults.

## Files

| File | Purpose |
|---|---|
| `compare_newton_paper.m` | Symmetric Newton CP comparison; internally implements RC-TLR node selection and a ChebTuck coefficient/RHOSVD comparator |
| `rc_tlr_1d_experiment.m` | Standalone RC-TLR experiment and results table |
| `rc_tlr_build.m` | General CP-input construction including embedded strong RRQR |
| `kernel_tensor.m` | Load original Newton data and form the requested full/long/short CP input, slice, or dense tensor |
| `cheb_tensor.m`, `my_chev2.m` | Chebyshev node sampling and spline evaluation |
| `cheb_weight_matrices.m`, `weighted_cheb_tensor.m`, `weighted_core.m` | Quadrature weighting |
| `selected_core.m`, `lagrange_cheb.m` | Unweighted selected core and barycentric Lagrange basis |
| `rc_htlr_cp_values.m`, `rc_htlr_relative_kernel_error.m` | CP block evaluation and error validation |
| `rc_htlr_1d_experiment.m` | Backward-compatible alias for the renamed entry point |

The grid storage reported is the selected/compressed core plus the three evaluated grid factors. The comparison driver reports median build times after a warm-up; the standalone driver measures construction of the selected core and separately retained interpolation and basis matrices. Neither timing includes reconstruction or evaluation of the reference.

## Attribution and verification

Benchmark data and the ChebTuck comparison methodology originate from Benner, Khoromskij, Khoromskaia and Sun, *A Hybrid Chebyshev--Tucker Tensor Format for Approximation of Multivariate Functions*, [arXiv:2503.01696](https://arxiv.org/abs/2503.01696). See the [original repository](https://github.com/bonans/ChebTuck) and [upstream MIT license notice](LICENSE-ChebTuck.txt). The standalone comparison is an implementation of the cited approach; it is not a claim that it is the unmodified upstream ChebTuck driver.

The MATLAB source files are checked in as supplied. Dependencies were inspected statically; MATLAB execution and verification against the manuscript tables are still pending. Input datasets and measured result tables are not included in this commit.
