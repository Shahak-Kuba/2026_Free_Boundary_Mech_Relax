# 2026_Free_Boundary_Mech_Relax

Simulation code for mechanical relaxation of a 1D cell population with a **free
boundary**, comparing the discrete cell-based model against its continuum limit.

The repository contains the Julia package `FreeBoundaryMechRelax` plus the scripts
that generate the paper figures. 

## The models

**Discrete.** `N` cells, each discretised into `m` springs (so `M = N·m` springs total), relaxing under overdamped Hookean or nonlinear spring mechanics. The left end is pinned at the origin; the right end is a free boundary at `x = L(t)`. Integrated with forward Euler at fixed step `δt`. For this project we ignore cell proliferation, death and embedment (flags are set by default to false) as that code was developed for alternate project.

**Continuum.** The corresponding free boundary PDE for the cell density, solved on the fixed computational domain `z ∈ [0,1]` with `x = L(t) z`, in two forms:

| Function | Boundary treatment |
|---|---|
| `rhs_Baker_in_ρ!` | No correction term at the free boundary (Baker et al. 2019) | | `rhs_with_correction!` | Includes the discrete correction term `\|C(t)\|` at the free boundary, evaluated by `L_series` |

`L_series(t, k, η, N, l0, a)` is the exact eigenfunction series for the correction
term; `compute_mF_q_term` is its leading-order short-/long-time asymptotic
approximation.

## Layout

```
src/
  FreeBoundaryMechRelax.jl              module definition
  Discrete/
    DataStructs.jl                      parameter and result types
    GeneralEquations.jl                  geometry helpers
    CellMechanics.jl                     spring force law F
    CellBehaviours.jl                    cell event probabilities, RNG seeding, saving callbacks
    EventCallbacks.jl                    proliferation / death / embedment callback
    ProblemSetup.jl                      initial conditions and ODEProblem construction
    PostSimulation.jl                    post-processing into SimResults_t
    FreeBoundary/
      FreeBoundaryODEProblem.jl          discrete RHS (FB_ODE!)
      FreeBoundarySimulation.jl          simulation driver
  Continuum/
    FreeBoundary/
      FB_1D_Solver.jl                    continuum RHS variants and parameters
      _superseded_rhs_variants.jl        archive, NOT included — see file header
scripts/
  figures_for_paper.jl                   reproduces the main paper figure
figures/                                 script output (git-ignored)
test/runtests.jl                         smoke tests
```

## Reproducing the figures

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia --project=. scripts/figures_for_paper.jl
```

This writes `figures/fig_for_paper_v1.pdf` (the two-panel main figure) and `figures/fig_left_zoomed.pdf` (the early-time inset). 

Run the tests with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

`Manifest.toml` is committed and carries the exact dependency versions used for
the published figures (Julia 1.10.5).

