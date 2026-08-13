using FreeBoundaryMechRelax
using OrdinaryDiffEq
using Test

const FBMR = FreeBoundaryMechRelax

@testset "FreeBoundaryMechRelax.jl" begin

    k_star, eta_star, a_star = 4.0, 1.0, 5.0
    L0, q₀, N_cells = 10.0, 1.0, 10
    T = 0.01   # short

    @testset "discrete free boundary simulation" begin
        Domain = FBMR.DomainProperties_t(N=N_cells, m=2, domain_type="1D")
        CellMech = FBMR.CellMechProperties_t(kₛ=k_star, η=eta_star, kf=0, a=a_star,
                                             restoring_force="hookean")
        SimTime = FBMR.SimTime_t(Tmax=T, δt=0.0001, event_δt=0.0001)
        FB_IC = FBMR.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = L0)
        ev = FBMR.CellEvent_t()

        sol = FBMR.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime,
                                          ev, ev, ev, ev, 1, 11)

        # left end stays pinned at the origin, and the free boundary starts at L0
        @test sol.u[1][1,1] == 0.0
        @test sol.u[end][1,1] == 0.0
        @test sol.u[1][1,end] ≈ L0 atol=1e-8

        # cells are longer than their resting length a*, so the tissue expands
        @test sol.u[end][1,end] > sol.u[1][1,end]

        # no events fire, so the cell count is conserved
        @test all(sol.CellCount .== N_cells)

        # node positions stay ordered
        @test all(diff(sol.u[end][1,:]) .> 0)
    end

    @testset "continuum free boundary solvers" begin
        k, η, a = k_star, eta_star, a_star
        α = k/η
        N = 101
        Δx = 1/(N-1)
        N_IC = Int(L0 * q₀)
        tspan = (0.0, T)

        Ffunc = ρ -> k * (1/ρ - a)
        Dfunc = ρ -> α / ρ^2
        Pfunc = ρ -> 0.0
        Afunc = ρ -> 0.0

        # without the discrete correction term (Baker et al. 2019)
        p_baker = FBMR.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed,
                                F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
        y0 = FBMR.make_initial_condition_FB(p_baker.N; U0fun = z -> q₀, L0=p_baker.L0)
        @test length(y0) == N + 1
        @test y0[end] == L0

        sol_baker = solve(ODEProblem(FBMR.rhs_Baker_in_ρ!, y0, tspan, p_baker),
                          Rodas5P(), saveat=[0.0, T])
        @test sol_baker.retcode == ReturnCode.Success
        @test sol_baker.u[end][end] > L0

        # with the discrete correction term
        p_corr = FBMR.FBParams_w_correction(α, k, a, η, Δx, L0, N, N_IC, q₀, :free,
                                            Ffunc, Dfunc, Pfunc, Afunc)
        sol_corr = solve(ODEProblem(FBMR.rhs_with_correction!, y0, tspan, p_corr),
                         Rodas5P(), saveat=[0.0, T])
        @test sol_corr.retcode == ReturnCode.Success
        @test sol_corr.u[end][end] > L0

        # in the m → ∞ limit the boundary node is frozen
        @test sol_corr.u[end][N] == sol_corr.u[1][N]
    end

    @testset "L_series" begin
        # |C(t)| decays monotonically and matches the closed form of its own sum
        vals = [FBMR.L_series(t, k_star, eta_star, N_cells, 1/q₀, a_star; P=300)
                for t in [0.1, 1.0, 10.0, 100.0]]
        @test all(diff(abs.(vals)) .< 0)

        # at long times a single mode dominates
        A_inf = (2 / (N_cells * eta_star)) * k_star * (1/q₀ - a_star)
        α = k_star / eta_star
        t = 500.0
        @test FBMR.L_series(t, k_star, eta_star, N_cells, 1/q₀, a_star; P=300) ≈
              A_inf * exp(-α * pi^2 * t / (4 * N_cells^2)) rtol=1e-6
    end
end
