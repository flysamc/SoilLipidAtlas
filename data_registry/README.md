# Data registry

Files too large for git do not live in this repository. Each row of
`registry.csv` records one such file:

| Column | Meaning |
|---|---|
| `file` | canonical name of the artifact |
| `step` | which methods step produced / uses it |
| `what` | one-line description (with key dimensions) |
| `size_bytes` | exact size |
| `sha256` | checksum of the released file |
| `current_location` | path used when the file was registered |
| `public_location` | where to obtain it: MassIVE, GNPS2 task IDs, or the Zenodo deposit |

To verify that a downloaded copy is exactly the file the analysis used:

```bash
sha256sum <file>   # must match the sha256 column
```

Public homes:

- **Raw spectra (.mzML):** MassIVE **MSV000102115**.
- **Per-batch FBMN results:** GNPS2 (task IDs in
  `methods/02_features/fbmn_batches_{POS,NEG}.csv`).
- **Processed large tables:** versioned, citable Zenodo deposit
  (DOI 10.5281/zenodo.20811187, reserved; populated from this registry).
