# Provenance

This merged source tree was prepared from the following author-supplied
artifacts. SHA-256 hashes identify the exact review inputs without embedding
machine-specific paths.

| Artifact | SHA-256 |
|---|---|
| `pStokesFEM.zip` | `cae294a05e084001748bb4b34d3ed886560cec474ed05d12f448fd0d3c86232f` |
| `pStokesFEMmixed.zip` | `99f69608123f309e73a5a70f8c2ad6881bf06b7d6e0c3c4027cc69e73bd26f81` |
| `PressureRate-21.pdf` | `c4bdc98f742dad653134e4f933f623a25f6d5c7c32c55a58afe0c5ca7c09d5ea` |
| `projections_code.zip` | `e1c950c187de52db00b309502137e195a3cb59f229e58b2111ee6c81ec335101` |
| `PressureRate-23.pdf` | `c09d8977ec046f71be642b6b79477960b4db5dc4b3eea524a4bf045306234eb6` |
| `PressureRate-36.pdf` | `6101daf6ff81ea88f111493940fad753a91ef5aac3105749ca04e2364e69c8a4` |

The original archives are not redistributed in this repository. Their two
implementations were reviewed against Definition 3.11 and the
numerical-experiment section of the supplied article. The principal behavioral
decisions retained or changed by the merge are documented in
`IMPLEMENTATION_NOTES.md`.

The numerical data included in `published_results/` were generated with the
merged code and are the data used for the tables and stability figure in the
associated article. The original pickle outputs from the supplied archives are
not redistributed.

For Section 6.2, `published_results/section_6_2_convergence/` contains all
64 level archives and eight EOC CSV files. For Section 6.3,
`published_results/section_6_3_projection_stability/` contains the stability
CSV and the generated strong/multiplier figure panels. The original plotting
script displayed only the strong-condition half of the figure; the packaged
plot entry point can generate both boundary-condition variants.
