"""
    FreeBoundaryMechRelax

Discrete and continuum models for mechanical relaxation of a 1D cell population
with a free boundary.

The package provides

* a **discrete** model — `N` cells, each discretised into `m` springs, relaxing
  under overdamped spring mechanics with a pinned left end and a free right end
  ([`FreeBoundarySimulation`](@ref)); and
* the corresponding **continuum** free boundary PDE, solved either without a
  correction term at the moving boundary (Baker et al. 2019,
  [`rhs_Baker_in_ρ!`](@ref)) or with the discrete correction term
  ([`rhs_with_correction!`](@ref), [`L_series`](@ref)).

This code is the subset of `MechCellTissueGrowth.jl` required to reproduce the
figures in the accompanying paper; see `scripts/figures_for_paper.jl`.
"""
module FreeBoundaryMechRelax

    # PACKAGES USED for solving equations
    using Base
    using OrdinaryDiffEq
    using DiffEqCallbacks
    using LinearAlgebra
    using Random
    using ElasticArrays
    using QuadGK
    using NLsolve

    # DEVELOPED SIMULATION CODE

    # shared discrete model code
    include("Discrete/DataStructs.jl")
    include("Discrete/GeneralEquations.jl")
    include("Discrete/CellMechanics.jl")
    include("Discrete/CellBehaviours.jl")
    include("Discrete/EventCallbacks.jl")
    include("Discrete/ProblemSetup.jl")
    include("Discrete/PostSimulation.jl")

    # discrete free boundary model
    include("Discrete/FreeBoundary/FreeBoundaryODEProblem.jl")
    include("Discrete/FreeBoundary/FreeBoundarySimulation.jl")

    # continuum free boundary model
    include("Continuum/FreeBoundary/FB_1D_Solver.jl")
    # NOTE: Continuum/FreeBoundary/_superseded_rhs_variants.jl is deliberately
    # NOT included — see the header of that file.

end
