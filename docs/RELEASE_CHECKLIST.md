# Zenodo release checklist

## Required author decisions

- [x] Use the BSD-3-Clause license.
- [x] Add the SPDX license identifier to `pyproject.toml`, `CITATION.cff`, and
      `.zenodo.json`.
- [x] Record the Zenodo DOI in `CITATION.cff` and the documentation.
- [x] Replace author initials with full names.
- [ ] Add ORCID identifiers to `CITATION.cff` and `.zenodo.json`, if desired.
- [x] Add `repository-code` to `CITATION.cff`.
- [ ] Add `related_identifiers` to `.zenodo.json` once the article DOI is known.
- [x] Confirm release version `1.0.0` / Git tag `v1.0.0`.

## Numerical verification

- [ ] Create the environment from `environment.yml`.
- [ ] Confirm that PETSc reports MUMPS availability.
- [ ] Run `pytest -q`; no FEniCS test may be skipped in the release log.
- [ ] Run the level-0 strong and multiplier smoke cases.
- [ ] Run `experiments/section_6_2_convergence/run.py`.
- [ ] Generate the manuscript table rows with
      `python scripts/generate_latex_tables.py --input results/section_6_2_convergence`.
- [ ] Run `experiments/section_6_3_projection_stability/run.py` for all 39
      meshes and generate the strong/multiplier figures.
- [ ] Repeat Section 6.3 on meshes 1, 8, 20, and 39 with quadrature degrees
      6, 10, and 14; document the maximum relative change.
- [ ] Confirm that `published_results/` contains the final data used in the
      article.
- [ ] Record compiler, PETSc, MUMPS, Python, NumPy, and platform versions.
- [ ] Decide whether the publication reproduces `--load-rule right` or also
      presents a time-averaged load experiment.

## Repository hygiene

- [ ] Remove generated caches and local results not intended for publication.
- [ ] Confirm that no pickle or machine-specific absolute path remains.
- [ ] Confirm that every command in `README.md` works from a clean checkout.
- [ ] Review `.zenodo.json` and `CITATION.cff` against Zenodo's preview.
- [ ] Set the release version consistently in `pyproject.toml`,
      `src/pstokes_fem/__init__.py`, `CITATION.cff`, and `CHANGELOG.md`.
- [ ] Create a GitHub release from tag `v1.0.0`.
- [x] Archive the exact Git tag on Zenodo and record the generated DOI.
