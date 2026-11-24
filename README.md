# RC-TLR
Quadrature-weighted interpolatory Tucker approximation using strong RRQR.

## Manuscript figures

[View the complete manuscript figure gallery](figures/README.md) (Figures 1–11, all 17 original image files, arranged by manuscript page and panel). The archived comparison figures are historical paper results; the executable MATLAB implementations remain RC-TLR-only.

## Reproduce all numerical experiments

Clone or download this repository, open its root directory in **MATLAB R2024b**, and run:

```matlab
summary = reproduce_all();
```

This executes every RC-TLR experiment currently included in this repository,
sequentially, using the bundled inputs. No ChebTuck implementation, Chebfun,
or Tensor Toolbox is needed. The biomolecular dataset files stored in
binary parts are restored and verified automatically (standard MATLAB JVM
required).

Start with the reduced real-data smoke run if you are using a new installation:

```matlab
summary = reproduce_all('quick');
```

The quick run verifies the execution paths on smaller cases; it is not a
reproduction of the full numerical tables. For a new output location,
including when changing the code or experiment settings, use an explicit path:

```matlab
summary = reproduce_all('full','C:\RC_TLR_results\run_1');
```

### Outputs

By default, results are written to `results/full/` (or `results/quick/`).

| Study | Full-run result files in its subfolder |
|---|---|
| Newton full and long-range tensors | `newton/newton_rc_tlr.csv`, `newton/newton_first_passing_p.csv`, individual grid tables, `grid_status.csv` |
| Standalone Newton CP-input RC-TLR | `newton_cp/newton_cp_rc_tlr.csv`, `newton_cp_results.mat` |
| Weighted/unweighted three-box RC-TLR | `three_box/threebox_rc_tlr.csv`, `threebox_weighting_comparison.csv`, generated error, storage and construction-time plots |
| Biomolecular (CP and full-grid spline) and lattice RC-TLR | `biomolecular/all_cases.csv`, study-specific tables, plots, `results.mat` |

`reproduction_summary.csv` records each stage's completion status, successful
row count, wall-clock elapsed time and result directory. The wall-clock time
includes the whole stage; it is separate from the experiment-specific
construction timings. A stage marked OK completed computation; that does not
mean every row satisfies the target error. Inspect each results table for
its accuracy checks. Failed stages or grids are recorded explicitly rather
than being represented as measured zero errors.

The full runner preserves the supplied numerical settings: the Newton
drivers use RRQR tolerance 5e-8, the three-box and biomolecular drivers
use 4e-8, and the target final relative error is 1e-7. The Newton sweeps
cover five grid sizes (70 total full/long-range configurations); the
standalone CP experiment uses its existing stopping criterion. The
three-box sweep has 14 weighted/unweighted configurations, and the
biomolecular suite has 55 configurations including the lattice. See
the module documentation for the exact error norms and storage definitions.
The full sweeps can take considerable time and memory; the quick profile
is recommended first, particularly on an 8 GB machine.

Stages run independently. Completed Newton grid tables are written as
each grid finishes; failures do not prevent subsequent stages from being
attempted. Biomolecular results use the module's existing checkpoint
mechanism. Re-running in an existing output directory reruns the other
stages; use a new directory for a clean run or after modifying the model.
No generated results are committed to Git.

**Verification:** the entry-point wiring and source interfaces were
inspected statically. MATLAB execution of the new runner remains to be
verified in a local MATLAB installation.

## Experiments

- [Coulomb / Newton kernel](experiments/newton/README.md): full and range-separated Newton tensors, RC-TLR approximation.
- [Effect of quadrature weighting on a three-box kernel](experiments/three_box/README.md): weighted and unweighted RC-TLR.
- [Biomolecular electrostatic potentials](experiments/biomolecular/README.md): standalone RC-TLR, including CP-input, full-grid spline and optional lattice studies.

See each module's README for its inputs, settings, and individual run instructions.
