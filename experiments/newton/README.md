# Coulomb / Newton kernel experiments

This module contains the supplied MATLAB implementations for the cell-averaged Coulomb (Newton) kernel on ([-1,1]^3), with both the full canonical tensor and its long-range component. It also contains a standalone RC-TLR CP-input experiment. Only RC-TLR is implemented; the strong RRQR tolerance follows the final manuscript.

## Requirements

- MATLAB (R2024b is the development environment).
- The Newton canonical factor data from [ChebTuck](https://github.com/bonans/ChebTuck/tree/2aa3729fe9991377bf1a6815ab3edd73423689dd/data/Newton_CP_-1_1).
- The Newton driver uses MATLAB built-in numerical routines; it does not call Chebfun or Tensor Toolbox.

## Data

Download the original `Newt_canon_<n>.mat` files from the pinned [ChebTuck dataset directory](https://github.com/bonans/ChebTuck/tree/2aa3729fe9991377bf1a6815ab3edd73423689dd/data/Newton_CP_-1_1). Each file contains the canonical factor matrix `U1r`. Keep these files under `experiments/newton/data/Newton_CP_-1_1/`, or pass the absolute path to your existing dataset directory.

Available grid sizes: 256, 512, 1024, 2048, 4096, 8192. The benchmark scales `U1r` by (1/h), where (h=2/n), and uses the cell-centred grid. The long-range part retains the first `ceil(R/2)-3` canonical columns, with (R) the full CP rank. These are cell-averaged kernel data, not point samples at the singularity.

## Quick start

From the repository root in MATLAB:

```matlab
newtonDir = fullfile(pwd,'experiments','newton');
addpath(newtonDir);
dataDir = fullfile(newtonDir,'data','Newton_CP_-1_1');

% Small RC-TLR run: both full and long-range parts, n=256, p=129.
[results, firstPass] = rc_tlr_newton_experiment(dataDir,256,129,129);
disp(results);
disp(firstPass);
```

To run the supplied full Newton sweep:

```matlab
[results, firstPass] = rc_tlr_newton_experiment(dataDir);
writetable(results,'newton_rc_tlr.csv');
writetable(firstPass,'newton_first_passing_p.csv');
```

For the standalone RC-TLR CP-input experiment:

```matlab
opts = struct('n',256,'pValues',129,'dataDir',dataDir, ...
              'warmups',0,'repeats',1,'qrTol',5e-8);
[results,model,diagnostics] = rc_tlr_1d_experiment(opts);
disp(results);
```

Increase node counts only after checking available memory. The general builder constructs a pairwise CP-factor product with `p^2` rows, whereas the symmetric Newton driver works in reduced orthonormal coordinates. Avoid dense `kernel_tensor(n)` calls for large grids; request `'cp'` or `'slice'` instead.

## Supplied numerical settings

| Setting | Symmetric Newton driver | Standalone RC-TLR |
|---|---|---|
| Grid sizes | 256, 512, 1024, 2048, 4096 | 256 |
| Full kernel node counts | 129, 257, 513, 1025, 2049, 4097, 8193, 16385 | Not the default input |
| Long-range node counts | 129, 257, 513, 1025, 2049, 4097 | 129, 257, 513, 1025 |
| RC-TLR RRQR tolerance | 5e-8 | 5e-8 (configurable) |
| Target relative error | 1e-7 | 1e-7 |
| Strong RRQR strength | 2 | 2 (configurable) |
| Validation | Relative maximum error, complete middle slice | Same by default; full-grid Frobenius error optional |

The RRQR tolerance and the target final approximation error are distinct parameters. Both RC-TLR Newton entry points use a default strong RRQR tolerance of 5e-8, as selected in the final manuscript. The target final error remains 1e-7.

## Files

| File | Purpose |
|---|---|
| `rc_tlr_newton_experiment.m` | RC-TLR-only symmetric full and long-range Newton sweeps |
| `rc_tlr_1d_experiment.m` | Standalone RC-TLR experiment and results table |
| `rc_tlr_build.m` | General CP-input construction including embedded strong RRQR |
| `kernel_tensor.m` | Load original Newton data and form the requested full/long/short CP input, slice, or dense tensor |
| `cheb_tensor.m`, `my_chev2.m` | Chebyshev node sampling and spline evaluation |
| `cheb_weight_matrices.m`, `weighted_cheb_tensor.m`, `weighted_core.m` | Quadrature weighting |
| `selected_core.m`, `lagrange_cheb.m` | Unweighted selected core and barycentric Lagrange basis |
| `rc_htlr_cp_values.m`, `rc_htlr_relative_kernel_error.m` | CP block evaluation and error validation |
| `rc_htlr_1d_experiment.m` | Backward-compatible alias for the renamed entry point |

The grid storage reported is the selected/compressed core plus the three evaluated grid factors. The Newton sweep reports core dimensions, Tucker rank, relative maximum error, core storage and total storage, without runtime columns. The separate standalone driver retains its original diagnostic build timing.

## Attribution and verification

Benchmark input data originate from Benner, Khoromskij, Khoromskaia and Sun, *A Hybrid Chebyshev--Tucker Tensor Format for Approximation of Multivariate Functions*, [arXiv:2503.01696](https://arxiv.org/abs/2503.01696). See the [original repository](https://github.com/bonans/ChebTuck) and [upstream MIT license notice](LICENSE-ChebTuck.txt). No ChebTuck approximation implementation is included.

The RC-TLR-only Newton driver retains the original RC-TLR construction and validation settings. Dependencies were inspected statically; MATLAB execution and verification against the manuscript tables are still pending. Input datasets and measured result tables are not bundled.
