# Publishing workflow

This repository is designed to be published as the reproducibility package for
the numerical experiments in Section 6 of the associated article.

## Recommended route

Use both GitHub and Zenodo:

1. Create a public GitHub repository, for example
   `pstokes-slip-pressure-code`.
2. Push exactly this directory as the repository root.
3. Replace the placeholder GitHub URL in `CITATION.cff`.
4. Replace placeholder DOI entries in `.zenodo.json` once known.
5. Confirm the BSD-3-Clause license metadata in GitHub and Zenodo.
6. Run the checks listed in `docs/RELEASE_CHECKLIST.md`.
7. Create a GitHub release, e.g. `v1.0.0`.
8. Archive that release on Zenodo and use the version-specific Zenodo DOI in
   the article.

GitHub provides a readable development repository; Zenodo provides the citable
and archived snapshot requested by many journals and referees.

## Suggested repository description

> Reproducibility code and published numerical data for finite-element experiments on
> pressure error estimates for the unsteady p-Stokes equations with slip
> boundary conditions.

## Suggested code-availability paragraph

Replace the placeholders after the Zenodo release has been minted:

```latex
\paragraph{Code availability.}
The legacy-FEniCS implementation and numerical data used for the numerical
experiments in Section~6 are archived on Zenodo
(version~v1.0.0, DOI: \url{https://doi.org/10.5281/zenodo.XXXXXXX}).
The development repository is available at
\url{https://github.com/OWNER/pstokes-slip-pressure-code}.
```

If the journal prefers a non-paragraph statement:

```latex
\noindent\textbf{Code availability.}
The code and numerical data used for the numerical experiments are available
from Zenodo at \url{https://doi.org/10.5281/zenodo.XXXXXXX}; the corresponding
development repository is
\url{https://github.com/OWNER/pstokes-slip-pressure-code}.
```

## License note

The code is distributed under the BSD-3-Clause license. This is a permissive
open-source license that allows reuse with attribution while retaining a
non-endorsement clause.
