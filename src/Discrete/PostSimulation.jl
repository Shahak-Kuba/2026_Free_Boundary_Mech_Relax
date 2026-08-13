"""
    PostCalcs1D(u, p, CellMech_at_t)

Perform post-calculation for 1D simulation data.

This function computes the cell densities for a given state `u` and parameters `p`.
Force, velocity, stress and curvature are returned as zeros: they are only
meaningful for the 2D evolving-interface problem.

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `p`: A tuple of parameters used in the calculations.
- `CellMech_at_t`: The spring mechanical properties in force at this time point.

# Returns
A tuple containing the sum of forces, normal velocity, density, stress, and curvature for each element in the state vector.

!!! note
    The original in `MechCellTissueGrowth.jl` also built shifted copies `uᵢ₊₁`,
    `uᵢ₋₁` and applied a periodic offset to them for the `"InvertedBellCurve"` and
    `"CosineSineWave"` domain types. Every use of those locals is commented out
    there, so they are dropped here; the returned values are unchanged.
"""
function PostCalcs1D(u, p, CellMech_at_t)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p

    kₘ = CellMech_at_t[1]
    aₘ = CellMech_at_t[2]
    kfₘ = CellMech_at_t[3]
    growth_dir = CellMech_at_t[4]
    η = CellMech.η

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    density = 1 ./ (Domain.m .* diff(u[1,:]))

    return ∑F, vₙ, density, ψ, Κ
end

"""
    postSimulation(sol, p, AllCellMech)

Perform post simulation calculations and return a comprehensive data structure with all relevant data.

# Arguments
- `sol`: The solution object from the simulation.
- `p`: Parameters used in the post calculations.
- `AllCellMech`: The saved spring mechanical properties at each save time.

# Returns
An instance of `SimResults_t` containing the calculated data.
"""
function postSimulation(sol, p, AllCellMech)

    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p

    c = size(sol.t, 1)
    s = min(c, size(AllCellMech,1))
    Area = Vector{Float64}(undef, s)
    Cell_Count = Vector{Float64}(undef, s)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Vector{Float64}}(undef, 0)
    DENSITY = Vector{Vector{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Vector{Float64}}(undef, 0)

    u = [Matrix((reshape(vec, 2, Int(length(vec)/2)))) for vec in sol.u]


    for ii in 1:s
        Area[ii] = Ω(u[ii]) # area calculation
        Cell_Count[ii] = (size(u[ii],2)-1)/Domain.m
        Fnet, nV, den, stre, kap = PostCalcs1D(u[ii], p, AllCellMech[ii])
        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, vec(den))
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(Domain.btype, sol.t[1:s], u[1:s], ∑F, DENSITY, vₙ, Area, ψ, Κ, Cell_Count)
end
