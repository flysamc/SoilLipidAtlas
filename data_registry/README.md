# Data registry

Big files don't live in git. Each row of `registry.csv` records one file: what
it is, which step uses it, its size, its SHA-256 checksum, where it currently
sits, and where it will be public (MassIVE / Zenodo).

To verify a downloaded copy is exactly the file the analysis used:

```bash
sha256sum <file>   # must match the sha256 column
```

Raw spectra (.mzML) are already public on MassIVE. Processed big tables will be
deposited on Zenodo as a versioned, citable package.
