# RC-TLR input datasets

All nine biomolecular CP datasets and the optional lattice input supplied
for this research project are included. The manifest records the original
byte counts, SHA-256 checksums, grid sizes, particle counts and CP ranks.

Three smaller files are stored directly. The seven larger files are stored
as ordered binary parts in `parts/`, because the upload interface could
not transfer the largest files as a single request. The parts are raw bytes,
not independent MAT files. From the module directory, run:

```matlab
restore_biomol_data();
```

This recreates and verifies the original files without changing any
existing file. `run_biomol_paper` automatically restores its requested
datasets before starting. Restored files are ignored by Git to avoid
committing redundant copies. An external folder with the original MAT files
can be passed as `dataRoot` instead.

For the lattice, the manifest reports the original CP rank 50; the supplied
long-range selection retains 28 terms. The biomolecular inputs are already
long-range CP tensors. No atomic coordinates or original preprocessing are
required by the approximation driver.
