# Changelog

## 1.0.0 - 2026-06-30

- First public reproducibility release for the numerical experiments in
  Section 6 of the associated article.
- Merge the strong and multiplier implementations around a shared residual.
- Preserve inhomogeneous normal-velocity and natural-traction data.
- Parameterize the pressure amplitude by the power-law exponent.
- Add optional temporal quadrature for the load and traction.
- Replace pickle output with compressed NumPy archives and JSON metadata.
- Add RT multiplier nullspace diagnostics and the published EOC tables.
- Document the lifted multiplier and provide reconstruction of the physical
  multiplier from Definition 3.11.
- Add command-line, packaging, citation, and Zenodo metadata.
- Add the Section 6.3 discrete-Leray projection experiment and its complete
  projection-stability figure workflow.
- Reuse assembled projector matrices and factorizations across basis vectors.
- Organize Section 6.2 and Section 6.3 as separate experiment folders with a
  shared implementation package.
- Add a LaTeX-table generator for the Section 6.2 EOC tables.
