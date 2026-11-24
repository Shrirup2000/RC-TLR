# RC-TLR: biomolecular electrostatic potentials

Standalone MATLAB implementation of quadrature-weighted, interpolatory
Tucker approximation using deterministic strong rank-revealing QR. The module
also retains the optional lattice experiment from the supplied numerical suite.
It has no baseline approximation implementation or baseline experiment driver.

## Run in MATLAB R2024b

From the repository root:

```matlab
addpath('experiments/biomolecular','-begin');
test_biomol_rc_tlr();  % synthetic tests, no external datasets required
out = run_biomol_paper();
```

Only base MATLAB is required. The included dataset files are restored
automatically from their binary parts (standard MATLAB JVM required for
SHA-256 checking). See [data/README.md](data/README.md).
To use an existing dataset folder instead, pass its path:

```matlab
out = run_biomol_paper('C:\path\to\Biomol\data');
```

Configure the strong RRQR tolerance and final error target directly:

```matlab
opts = struct('qrTol',4e-8,'target',1e-7);
out = run_biomol_paper([],opts);
```

The default plan uses 55 distinct RC-TLR cases after deduplication, with
first-kind Chebyshev nodes (polynomial degree at most `p-1`).

| Study | Configuration | CSV |
|---|---|---|
| Grid refinement, CP input | N=500; n=256,512,1024,2048; p=129 | Table3.csv |
| Full-grid spline input | n=256; N=500; p=129 | Algorithm32.csv |
| Node and rank sweep | n=2048; N=100:100:600; p=9,17,33,65,129,257,513,1025 | Table4.csv, Ranks_RC_TLR.csv |
| Storage scaling | n=2048; N=100:100:600; p=129 | Table5.csv |
| Interpolation choice within RC-TLR | n=2048; N=500; p=129; spline, linear, nearest | Interpolation.csv |
| Optional lattice | n=960; N=2304; p=129 | Lattice.csv |

To run only the four biomolecular grid sizes, for example:

```matlab
opts = struct('qrTol',4e-8,'target',1e-7, ...
    'runRankStudy',false,'runInterpolation',false, ...
    'includeAlgorithm32',false,'runLattice',false,'plots',false);
out = run_biomol_paper([],opts);
```

To change the node counts of the rank study, set `pValues`; to change the
main node count, set `p`. The full-grid experiment is restricted to
`n<=256` and `p<=129` to limit memory, and larger CP sweeps can be slow.
`workspaceMB` sizes the block QR workspace (default 128 MiB); it is not a
bound on overall MATLAB process memory.

## Implementation

- `rc_tlr_biomol_experiment.m`: CP-input RC-TLR, reference-slice validation,
  full-grid route, outputs and checkpointing.
- `rc_tlr_functional.m`: direct approximation from a vectorized function.
- `rc_tlr_from_grid.m`: full-grid cubic spline followed by functional RC-TLR.
- `rc_tlr_functional_factors.m`: evaluate the compressed continuous model.
- `private/`: weighted sampling, strong RRQR, unweighted selected core,
  Lagrange basis and CP contractions. The block CP reduction is orthogonal,
  without a Gram-matrix rank approximation. The p=129 biomolecular branch
  retains the supplied dense weighted-sample route.
- `biomol_experiment_plan.m`: all defaults and study switches.

The biomolecular inputs provide CP weights `LRed.LAM1C/Hunif^3` and factors
`LRed.CU1`, `CU2` and `CU3`. Source coordinates use cell centers on
[-1,1], not endpoint-inclusive grids. The lattice loader retains the supplied
scaling, normalization and removal of eleven columns from each signed
component (50 original terms, 28 long-range terms).

## Accuracy and storage

Each experiment validates the **complete middle slice** at MATLAB index
`floor(n/2)`; it does not assert a full-volume error. If A is the original
CP middle slice, P the uncompressed Chebyshev interpolant and B the RC-TLR
approximation, the exported metrics are:

| Metric | Definition |
|---|---|
| AbsoluteError | max(abs(B(:)-A(:))) |
| RelativeError | AbsoluteError / max(abs(A(:))) |
| RelativeFrobeniusError | norm(B(:)-A(:)) / norm(A(:)) |
| CompressionError | max(abs(B(:)-P(:))) |
| InterpolationError | max(abs(P(:)-A(:))) |
| MaxQRResidual, MaxRho | Maximum weighted mode residual and strength statistic |
| CoreStorage | Dimensions rX × rY × rZ |
| MaxRank | max(rX,rY,rZ) |
| ContinuousStorage | rX*rY*rZ + p*(rX+rY+rZ) |
| GridStorage | rX*rY*rZ + n*(rX+rY+rZ) |
| FactoredGridStorage | ContinuousStorage + 3*n*p (CP biomolecular route) |

All storage counts are numbers of scalars, not peak RAM.
The strong RRQR tolerance is an internal criterion, not a guarantee that
the final slice error meets the target: `RelativeToleranceMet` and
`AbsoluteToleranceMet` record the observed checks.

Construction-only timing is retained as a diagnostic from the supplied
driver. For biomolecular CP input, the single timed construction includes
the core, extension maps and grid Lagrange bases. Other routes time continuous
construction. Factor products and reconstruction are separately recorded;
reference evaluation, diagnostics, loading and plotting are excluded from
`TotalTime`. The reported time is not a full end-to-end runtime.

The default results directory is `results/` beside this README, containing
`all_cases.csv`, study tables, plots, and `results.mat`. Failed cases include
error reports and are not measurements. Re-running resumes matching
successful cases. Change `outDir` or set `resume=false` after modifying settings
or source files. The repository includes no fabricated benchmark measurements.

See [VALIDATION.md](VALIDATION.md) for the status of checks.
