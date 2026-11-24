# Newton CP input data

The following original MAT files are included, with bytes verified against [the pinned ChebTuck repository](https://github.com/bonans/ChebTuck/tree/2aa3729fe9991377bf1a6815ab3edd73423689dd/data/Newton_CP_-1_1):

- `Newt_canon_256.mat`
- `Newt_canon_512.mat`
- `Newt_canon_1024.mat`
- `Newt_canon_2048.mat`
- `Newt_canon_4096.mat`
- `Newt_canon_8192.mat` (optional; not used by the default Newton sweep)

The original benchmark data are bundled here for reproducibility. Alternatively, pass their actual directory to the Newton experiment function or set `opts.dataDir` for the standalone experiment. See the [experiment instructions](../../README.md).
