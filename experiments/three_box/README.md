# Effect of quadrature weighting on a three-box kernel

This experiment isolates quadrature weighting within RC-TLR. Both variants apply strong RRQR to grouped mode unfoldings, using weighted and unweighted samples respectively. Both extract their interpolatory core from the unweighted sample tensor.

The kernel is \(K(x,y,z)=(\|x-y\|_2^2+\|y-z\|_2^2+\|z-x\|_2^2)^{-1/2}\). We use uniform 20 by 20 planar grids.

| Box | Domain |
|---|---|
| X | [-3,-1] x [0,2] |
| Y | [1,3] x [0,2] |
| Z | [-1,1] x [4,6] |

| Parameter | Default |
|---|---|
| Chebyshev nodes per coordinate | 2:8 (first kind) |
| Grouped mode size | p^2 |
| Strong RRQR strength | 2 |
| Strong RRQR tolerance | 4e-8 |
| Target relative Frobenius error | 1e-7 |
| Validation | Complete Cartesian product, evaluated in blocks |

## Run

From the repository root in MATLAB R2024b:

```matlab
addpath('experiments/three_box');
[results, comparison] = rc_tlr_threebox_experiment();
disp(results);
writetable(results,'threebox_rc_tlr.csv');
```

The driver contains all required local helpers. It produces relative error, ranks, compact storage, and construction time for weighted and unweighted RC-TLR. Construction timing is retained here because it is an explicit figure in this experiment; the Newton module does not record timing. For a Cartesian grid, storage counts the unweighted core, three interpolation maps, and six univariate basis matrices: rX*rY*rZ + p^2*(rX+rY+rZ) + 6*n*p. Reconstruction and reference evaluation are not included in construction time.

## Manuscript figures (original supplied plots)

The three supplied figures are preserved unchanged, including the previously measured ChebTuck baseline. The new driver generates only the weighted/unweighted RC-TLR curves; it does not execute or reproduce the ChebTuck baseline. The original figures are historical experimental outputs, not newly measured results from this refactored driver.

| Figure | File |
|---|---|
| Relative Frobenius error | [error.png](figures/error.png) |
| Scalars stored | [storage.png](figures/storage.png) |
| Construction time | [build.png](figures/build.png) |

The reference ChebTuck methodology is described in Benner, Khoromskij, Khoromskaia and Sun, [A Hybrid Chebyshev–Tucker Tensor Format for Approximation of Multivariate Functions](https://arxiv.org/abs/2503.01696); see the [original repository](https://github.com/bonans/ChebTuck). Its approximation implementation is not included in this module.

The RC-TLR construction and error routines were retained from the supplied experiment, with associated grouped-coordinate helpers recovered from the earlier experiment files. The refactored code has been inspected statically; MATLAB execution and numerical verification against the manuscript remain pending.
