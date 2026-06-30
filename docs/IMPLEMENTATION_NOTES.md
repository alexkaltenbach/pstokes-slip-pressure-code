# Implementation notes

## Changes relative to the supplied archives

- The two complete code copies are replaced by a shared package and two short
  boundary-condition functions.
- `PlotDG` and plotting imports are removed.
- Pressure amplitude is selected from the power-law exponent.
- Absolute and relative nonlinear tolerances are consistently set to
  `1e-10` and `1e-8`.
- The RT interior extension uses ordinary `dS`. In the supplied code, a
  marked measure was created but used without `(1)`, which also means
  `everywhere` in UFL; the expensive marker construction therefore had no
  effect.
- The fourth mixed component is consistently named `lifted_multiplier`.
  `boundary.article_multiplier` reconstructs the physical multiplier from
  Definition 3.11 by adding back the known normal stress and sign change.
- The default zero stress is represented as a tensor whenever one is needed.
- `raise NotImplemented` is replaced by explicit `ValueError` exceptions.
- Modules have no computation side effects on import.
- Raw results use `.npz` plus JSON metadata rather than pickle.
- Filenames contain no colons and are portable across operating systems.
- Compiler optimization defaults to portable `-O3`; `-ffast-math` and
  `-march=native` are not enabled by the packaged code defaults.
- The Section 6.3 projector matrices and direct factorizations are reused for
  all nodal basis vectors on a fixed mesh. The prototype rebuilt and refactored
  the same systems for every basis vector.
- Section 6.3 output is one checkpointed CSV instead of several pickle files;
  one traversal evaluates both values of `p`, both projectors, and all chosen
  quadrature degrees.

## Runtime verification status

The merged source is statically syntax checked in an environment without
DOLFIN. End-to-end numerical verification is performed in a FEniCS 2019.1
environment. The numerical data used in the associated article are included in
`published_results/`.

The projection indicator uses quadrature degree 6 by default, matching the
prototype. Non-integer powers are not polynomial integrands, so the Section
6.3 entry point also accepts several quadrature degrees in one run. A small
mesh subset at degrees 6, 10, and 14 is part of the recommended release check.


## Deliberately retained behavior

- inhomogeneous normal velocity and natural traction;
- right-endpoint load evaluation as the default reproduction mode;
- continuation paths used by the original implementations;
- Taylor-Hood plus global pressure mean multiplier;
- RT facet restriction in the weak-imposition variant;
- the intended `L-infinity-in-time` velocity error.
