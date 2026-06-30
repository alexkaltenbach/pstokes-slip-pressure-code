# Section 6.2: manufactured-solution convergence experiment

This folder contains the entry points for the convergence tables. The finite
element method itself lives in `src/pstokes_fem/` and is shared by both normal
boundary-condition implementations.

Run the complete parameter sweep with

```bash
python experiments/section_6_2_convergence/run.py \
  --output results/section_6_2_convergence
```

Level 7 uses 512 time steps and a uniformly refined Taylor--Hood mesh, so the
complete sweep is expensive. Existing result files are retained unless
`--overwrite` is passed.

The script writes one `*_eoc.csv` file for each parameter combination. The
published EOC tables shipped with the archive are in
`published_results/section_6_2_convergence`.

Generate the LaTeX rows used in Tables 1 and 2 from the shipped data with

```bash
python scripts/generate_latex_tables.py
```

or from newly generated data with

```bash
python scripts/generate_latex_tables.py \
  --input results/section_6_2_convergence
```
