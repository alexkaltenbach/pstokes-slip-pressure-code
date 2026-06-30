# Section 6.3: projection-stability indicator

This experiment evaluates

$$
c_{\mathrm{stab}}^i(\mathcal J_h)
=\max_j
\frac{\|\mathcal J_h P_{V_h}\varphi_j\|_{L^r}}
     {\|P_{V_h}\varphi_j\|_{L^r}},
\qquad
\mathcal J_h\in\{\mathcal P_h,\mathcal P_h^\perp\},
$$

for the nodal basis of the unconstrained vector-valued P2 space. It is the
basis-function indicator used in the article; it is not a proof of a uniform
operator-norm bound.

The default run covers the 39 meshes, both values of `p`, both projectors, and
both implementations of the normal condition:

```bash
python experiments/section_6_3_projection_stability/run.py
python experiments/section_6_3_projection_stability/plot.py \
  --input results/section_6_3_projection_stability/stability.csv \
  --enforcement strong \
  --output results/section_6_3_projection_stability/stability_strong.pdf
python experiments/section_6_3_projection_stability/plot.py \
  --input results/section_6_3_projection_stability/stability.csv \
  --enforcement multiplier \
  --output results/section_6_3_projection_stability/stability_multiplier.pdf
```

The run is serial and computationally expensive. Results are checkpointed
after each mesh/enforcement pair, so an interrupted run can simply be resumed.
Matrices and direct factorizations are reused while traversing the basis.

The published stability CSV and figure outputs shipped with the archive are in
`published_results/section_6_3_projection_stability`. The plotting script uses
that CSV by default, so the published figure panels can be regenerated with

```bash
python experiments/section_6_3_projection_stability/plot.py \
  --enforcement strong \
  --output results/section_6_3_projection_stability/stability_strong.pdf
python experiments/section_6_3_projection_stability/plot.py \
  --enforcement multiplier \
  --output results/section_6_3_projection_stability/stability_multiplier.pdf
```

The article used quadrature degree 6. Since `p`, `p'` can be non-integers, a
small sensitivity run is recommended:

```bash
python experiments/section_6_3_projection_stability/run.py \
  --mesh-indices 1,8,20,39 \
  --quadrature-degrees 6,10,14 \
  --output results/section_6_3_projection_stability/quadrature_check.csv
```

Use `--overwrite` to replace existing blocks. The CSV format replaces the
prototype's pickle output and records the mesh size, exponent, projector,
normal-condition implementation, and quadrature degree explicitly.
