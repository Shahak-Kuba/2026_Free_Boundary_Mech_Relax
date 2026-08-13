"""
    generate_homogeneous_population_FB(C1::CellMechProperties_t, N, m)

Expand a single cell type `C1` into a homogeneous population of `N` cells, each
discretised into `m` springs, for the free boundary problem. The `+ 1` accounts
for the extra node at the free boundary.

Spring stiffness is scaled by `m`, while the resting length, tissue production
rate and drag coefficient are scaled by `1/m`, so that the cell-level mechanics
are independent of `m`.
"""
function generate_homogeneous_population_FB(C1::CellMechProperties_t, N, m)
    M = (N * m) + 1
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{String}([])

    for i in 1:M
        ks[i] = C1.kₛ * m
        kf[i] = C1.kf / m
        a[i] = C1.a / m
        push!(growth_dir, C1.growth_dir)
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

# For free boundary Problem

"""
    generate_discrete_IC_from_density_profile(Func, Func_Derivative, M, L0, eps=0.1)

Place `M` spring boundaries on `[0, L0]` such that the resulting discrete density
matches the continuum density profile `Func` (with derivative `Func_Derivative`).
The interior node positions are found by solving the nonlinear system that
equates the discrete density gradient to `λ * Func'(x)`, where `λ` normalises
`Func` to `M` springs.

!!! note
    `MechCellTissueGrowth.jl` also defined a `(Func, Func_Derivative, N, m, L0,
    eps=0.1)` variant that computed `M = N*m` internally. Because both signatures
    are five untyped positional arguments, that method was overwritten by this one
    and was never reachable; it is not ported. (Keeping both also makes the module
    fail precompilation.) Pass `M = N*m` yourself.
"""
function generate_discrete_IC_from_density_profile(Func, Func_Derivative, M, L0, eps=0.1)
    # Initial condition setup
    integral, error = quadgk(x -> Func(x), 0, L0)
    lambda = M/integral

    function x_interior!(F,x)
        for i in 1:length(F)
            if i == 1
                F[i] = ((1/(x[i]-0)) + (1/(x[i]-x[i+1]))) / ((0-x[i+1])/2) - lambda*Func_Derivative(x[i])
            elseif i == length(F)
                F[i] = ((1/(x[i]-x[i-1]))+(1/(x[i]-L0))) / ((x[i-1]-L0)/2) - lambda*Func_Derivative(x[i])
            else
                F[i] = ((1/(x[i]-x[i-1]))+(1/(x[i]-x[i+1]))) / ((x[i-1]-x[i+1])/2) - lambda*Func_Derivative(x[i])
            end
        end
    end

    x_guess = collect(range(eps, stop=L0-eps, length=M-1))
    sol = nlsolve(x_interior!, x_guess)
    x_interior_sol = sol.zero
    initial_condition_discrete = vcat(0.0, x_interior_sol, L0)
    return initial_condition_discrete
end

"""
    SetupFBODEproblem(FB_IC, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

Build the `ODEProblem` for the free boundary discrete model from an initial
density profile `FB_IC`.
"""
function SetupFBODEproblem(FB_IC, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    q0 = FB_IC.q0
    q0_der = FB_IC.q0_der
    L0 = FB_IC.L0

    x0 = generate_discrete_IC_from_density_profile(q0, q0_der, M, L0, 0.001)
    u0 = ElasticMatrix([x0'; zeros(1, length(x0))])
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(FB_ODE!,u0,tspan,p), p
end

"""
    SetupFBODEproblem_given_IC(IC, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

Build the `ODEProblem` for the free boundary discrete model from an explicit
vector of initial node positions `IC`.
"""
function SetupFBODEproblem_given_IC(IC::Vector{Float64}, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

    u0 = ElasticMatrix(zeros(2,M+1))
    u0[1,2:end] = IC;
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(FB_ODE!,u0,tspan,p), p
end
