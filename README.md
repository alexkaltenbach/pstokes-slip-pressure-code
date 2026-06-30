# p-Stokes slip pressure experiments

Reproducibility package for the numerical experiments in Section 6 of

> Alex Kaltenbach and Jörn Wichmann,
> *A Priori Error Analysis for the p-Stokes Equations with Slip Boundary
> Conditions: A Discrete Leray Projection Framework*, 2026.

The repository is intended to accompany the article as an archived software
snapshot.  The suggested repository/archive name is
`pstokes-slip-pressure-code`; the Python import package remains
`pstokes_fem`.

The repository contains two discretizations of the impermeability condition:

- `strong`: the normal velocity is imposed component-wise on the
  axis-aligned unit square;
- `multiplier`: the normal velocity is imposed weakly with the normal trace of
  a facet-restricted Raviart-Thomas element. The stored unknown is a lifted
  multiplier; the physical multiplier from the article can be reconstructed
  with `pstokes_fem.boundary.article_multiplier`.

Both variants share the same p-Stokes operator and discrete Leray projectors.
The experiments from Sections 6.2 and 6.3 have separate entry points under
`experiments/`; mathematical implementation code is shared in `src/pstokes_fem/`.

## Mathematical scope

The manufactured experiment uses inhomogeneous boundary data

$$
u\cdot n=g_n, \qquad (S(Du)n)_\tau=g_\tau,
$$

because this permits explicit solutions with controlled fractional regularity.
See [`docs/NUMERICAL_METHOD.md`](docs/NUMERICAL_METHOD.md) for the precise
discrete forms and the multiplier lifting.

The default `right` load rule reproduces the supplied research code: volume
force and natural traction are evaluated at the right time endpoint. Optional
`midpoint`, `gauss2`, and `gauss3` rules approximate their interval averages
without symbolic time integration.

## Environment

The code targets FEniCS/DOLFIN 2019.1.0. The PyPI package named `fenics` does
not include DOLFIN and is therefore insufficient. A Conda environment is
provided:

```bash
conda env create -f environment.yml
conda activate pstokes-slip-pressure-code
python -m pip install -e ".[test]"
```

If a compatible FEniCS 2019.1 environment already exists, it can be reused:

```bash
conda activate fenicsproject
python -c "import dolfin; print(dolfin.__version__)"
python -m pip install -e ".[test]"
```

The PETSc installation must provide MUMPS, matching the solver used for the
reported experiments.

## Experiment 1: Section 6.2

Run a small strong-boundary experiment:

```bash
pstokes-run \
  --enforcement strong \
  --p 1.5 \
  --alpha 0.5 \
  --levels 0:2 \
  --output results
```

Run the multiplier variant:

```bash
pstokes-run \
  --enforcement multiplier \
  --p 2.5 \
  --alpha 1.0 \
  --levels 0:2 \
  --output results
```

Reproduce the complete convergence sweep from the paper:

```bash
python experiments/section_6_2_convergence/run.py
```

This is computationally expensive: level 7 contains 512 time steps and a
uniformly refined Taylor-Hood mesh.

Generate LaTeX rows for Tables 1 and 2 from a completed sweep:

```bash
python scripts/generate_latex_tables.py \
  --input results/section_6_2_convergence
```

Generate the LaTeX rows from the published data shipped with this archive:

```bash
python scripts/generate_latex_tables.py
```

See [`experiments/section_6_2_convergence/README.md`](experiments/section_6_2_convergence/README.md)
for the generated files.

## Experiment 2: Section 6.3

Evaluate the basis-function indicator for the discrete Leray projector and its
complement, for both treatments of the normal condition:

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

The projector systems are assembled and factorized once per mesh and then
reused for every basis vector. The run checkpoints after each mesh and normal
condition. See
[`experiments/section_6_3_projection_stability/README.md`](experiments/section_6_3_projection_stability/README.md)
for quadrature-sensitivity checks and computational details.

The same plotting script can be used directly with the published stability
data shipped in this archive:

```bash
python experiments/section_6_3_projection_stability/plot.py \
  --enforcement strong \
  --output results/section_6_3_projection_stability/stability_strong.pdf
python experiments/section_6_3_projection_stability/plot.py \
  --enforcement multiplier \
  --output results/section_6_3_projection_stability/stability_multiplier.pdf
```

## Results

Each run writes a compressed NumPy archive (`.npz`) containing the raw
per-time-step error integrals and a JSON metadata record. Unlike pickle, this
format does not execute code when loaded. The EOC summary is written as CSV.

The numerical data used for the tables and stability figure in the associated
article are included in [`published_results`](published_results/). New runs
write to [`results`](results/), which is ignored by Git.

## Tests

```bash
pytest -q
```

The pure-Python tests run without FEniCS. FEniCS-dependent tests are skipped
when DOLFIN is unavailable. In particular,
`tests/test_multiplier_nullspace.py` checks the rank of the RT trace mass
matrix and confirms that the original marked measure, used without a marker
id, is equivalent to ordinary `dS`. `tests/test_projections.py` checks the
discrete Leray decomposition and its L2 orthogonality on a small mesh.

## Repository layout

```text
src/pstokes_fem/                         reusable implementation
experiments/section_6_2_convergence/     convergence-table entry points
experiments/section_6_3_projection_stability/  projector experiment and plot
tests/                                   unit and diagnostic tests
published_results/                       numerical data used in the article
docs/                                    numerical method and release notes
results/                                 generated outputs (ignored by Git)
```

## Archiving and citation

The recommended publication workflow is:

1. create a public GitHub repository from this directory;
2. create a GitHub release, e.g. `v1.0.0`;
3. archive that release on Zenodo;
4. cite the versioned Zenodo DOI in the article.

See [`docs/PUBLISHING.md`](docs/PUBLISHING.md) for a step-by-step checklist
and a suggested code-availability paragraph.

The code is distributed under the BSD-3-Clause license; see [`LICENSE`](LICENSE).

The archived release is available on Zenodo:
[`10.5281/zenodo.21071175`](https://doi.org/10.5281/zenodo.21071175).

## Remaining release metadata

The remaining release metadata depends on optional author/article identifiers:

1. add ORCID identifiers to `CITATION.cff` and `.zenodo.json`, if desired;
2. add article DOI metadata to `.zenodo.json` once the article DOI is known;
3. confirm whether a metadata-only follow-up release is desired.

See [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md).
