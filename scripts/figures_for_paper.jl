import FreeBoundaryMechRelax as FBMR
using OrdinaryDiffEq
using QuadGK
using LinearAlgebra
using CairoMakie
using Roots, SpecialFunctions

# Figures are written next to the repository root, independent of the current
# working directory.
const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures"))
mkpath(FIGDIR)

### Colors for plotting
const BLACK          = "#000000"
const ORANGE         = "#E69F00"
const SKY_BLUE       = "#56B4E9"
const BLUISH_GREEN   = "#009E73"
const YELLOW         = "#F0E442"
const BLUE           = "#0072B2"
const VERMILLION     = "#D55E00"
const REDDISH_PURPLE = "#CC79A7"
###



N_cells      = 10            # the integer N in theta_p formula
k_star = 4          # k*
eta_star = 1          # eta*  (avoiding name collision with `eta`)
ell0   = 1.0          # l(0)
a_star = 5          # a*
ell0_minus_astar = ell0 - a_star
α  = k_star / eta_star

T = 50

# setting domain size
L0 = 10.0;
q₀ = 1.0;

# discrete model Params
CellMech = FBMR.CellMechProperties_t(kₛ=k_star, η=eta_star, kf = 0, a = a_star, restoring_force="hookean")
Prolif = FBMR.CellEvent_t()
Death = FBMR.CellEvent_t()
Embed = FBMR.CellEvent_t()
ProlifEmbed = FBMR.CellEvent_t()
SimTime = FBMR.SimTime_t(Tmax=T, δt=0.0001, event_δt=0.0001)
FB_IC = FBMR.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = L0)
NumSaveTimePoints = 1001

m_vals = [1,2,5,20]

# output from discrete simulations
all_disc_solutions = []
all_disc_solutions_boundary = []
t_disc = []

for m_value in m_vals
    # discrete model Params
    Domain = FBMR.DomainProperties_t(N=N_cells, m = m_value, domain_type="1D")
    # run discrete simulation
    sol_disc = FBMR.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    push!(all_disc_solutions, sol_disc)
    # data from discrete simulation
    solution_Boundary = zeros(length(sol_disc.t))
    for i in eachindex(sol_disc.t)
        solution_Boundary[i] = sol_disc.u[i][1,end]
    end
    push!(all_disc_solutions_boundary, solution_Boundary)
    push!(t_disc, sol_disc.t)
end


# Solving PDE
# Parameters
k = CellMech.kₛ; η = CellMech.η; α = k/η; a = CellMech.a; L0 = L0; N = 1001; Δx = 1/(N-1);
N_IC = Int(L0 * q₀); rBC=:free;
tspan = (0.0, SimTime.Tmax);

# Force and Diffusivity functions (Nonlinear spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2
# Proliferation and apoptosis functions
Pfunc = ρ -> 0.0
Afunc = ρ -> 0.0

q0_func = z -> q₀

# solving Baker et al. 2019 continuum model 
p = FBMR.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
y0 = FBMR.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)

# solve PDE
prob = ODEProblem(FBMR.rhs_Baker_in_ρ!, y0, tspan, p)
sol_Baker = solve(prob, Rodas5P(), saveat=vcat([0.0:0.1:tspan[2]]...))

# solving new continuum model with correction term
p = FBMR.FBParams_w_correction(α, k, a, η, Δx, L0, N, N_IC, q₀, rBC, Ffunc, Dfunc, Pfunc, Afunc)
y0 = FBMR.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(FBMR.rhs_with_correction!, y0, tspan, p)
sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.1:tspan[2]]...))

# plotting L(t) for discrete and continuum models
f = Figure(size=(1800,800));
ax = Axis(f[1,1], aspect=1, xlabel=L"$t$", ylabel=L"$L(t)$", limits=(-1,51,9,51));
L_cont = [sol_Baker.u[i][end] for i in eachindex(sol_Baker.t)]
lines!(ax, sol_Baker.t, L_cont, label=L"\text{Baker et al. 2019}", linewidth=4, color = BLUISH_GREEN);

# continuum solution for L(t) with correction term
L_cont = [sol.u[i][end] for i in eachindex(sol.t)]
lines!(ax, sol.t, L_cont, label=L"\text{Continuum limit}", linewidth=4, color = BLACK);
# discrete solution for L(t)
clrs = [VERMILLION, SKY_BLUE, REDDISH_PURPLE, YELLOW]
for ii in eachindex(all_disc_solutions)
    scatter!(ax, all_disc_solutions[ii].t[1:100:end], all_disc_solutions_boundary[ii][1:100:end], label=L"m=%$(m_vals[ii]),\; M = %$(N_cells*m_vals[ii])", markersize=15, color=clrs[ii], strokecolor=BLACK, strokewidth=0.8);
    scatter!(ax, all_disc_solutions[ii].t[end], all_disc_solutions_boundary[ii][end], markersize=15, color=clrs[ii], strokecolor=BLACK, strokewidth=0.8);
end
axislegend(ax, position=:rb, framevisible=false, patchsize=(40, 10))
display(f)


# subfigure for plotting convergence of discrete solution to continuum solution at early times
f2 = Figure(size=(250,250));
ax = Axis(f2[1,1], aspect=1, limits=(9,11,35,40), xticks=[9, 10, 11], yticks=[35, 40]);
L_cont = [sol_Baker.u[i][end] for i in eachindex(sol_Baker.t)]
lines!(ax, sol_Baker.t, L_cont, label=L"\text{Baker et al. 2019}", linewidth=5, color = BLUISH_GREEN);
# continuum solution for L(t) with correction term
L_cont = [sol.u[i][end] for i in eachindex(sol.t)]
lines!(ax, sol.t, L_cont, label=L"\text{Continuum limit}", linewidth=5, color = BLACK);
# discrete solution for L(t)
clrs = [VERMILLION, SKY_BLUE, REDDISH_PURPLE, YELLOW]
for ii in eachindex(all_disc_solutions)
    scatter!(ax, all_disc_solutions[ii].t[1:100:end], all_disc_solutions_boundary[ii][1:100:end], label=L"m=%$(m_vals[ii]), M = %$(N_cells*m_vals[ii])", markersize=15, color=clrs[ii], strokecolor=BLACK, strokewidth=0.8);
    scatter!(ax, all_disc_solutions[ii].t[end], all_disc_solutions_boundary[ii][end], markersize=15, color=clrs[ii], strokecolor=BLACK, strokewidth=0.8);
end
display(f2)
save(joinpath(FIGDIR, "fig_left_zoomed.pdf"), f2)

# -----------------------------------------------------------------------------
# Parameters (chosen arbitrarily; you can change these)
# -----------------------------------------------------------------------------
N      = 10            # the integer N in theta_p formula
k_star = 4          # k*
eta_star = 1          # eta*  (avoiding name collision with `eta`)
ell0   = 1.0          # l(0)
a_star = 5          # a*
ell0_minus_astar = ell0 - a_star
α  = k_star / eta_star

T = 50

# -----------------------------------------------------------------------------
# Running discrete simulation for m = 20
# -----------------------------------------------------------------------------

m = 20
q₀ = 1.0;
L0 = 10.0;
Domain = FBMR.DomainProperties_t(N=N_cells, m = m, domain_type = "1D")
CellMech = FBMR.CellMechProperties_t(kₛ=k_star, kf = 0, a = a_star, η = eta_star, restoring_force="hookean")
Prolif = FBMR.CellEvent_t()
Death = FBMR.CellEvent_t()
Embed = FBMR.CellEvent_t()
ProlifEmbed = FBMR.CellEvent_t()
SimTime = FBMR.SimTime_t(Tmax=T, δt=0.0001, event_δt=0.0001)
FB_IC = FBMR.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = L0)
Seed = 1
NumSaveTimePoints = 50001
sol = FBMR.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints);

idx_for_plots = [2, 11, 101, 1001, 10001, 20001, 30001, 40001, 50001]

first_term_disc = -[Domain.m / CellMech.η * CellMech.kₛ * (1/sol.Density[ii][end] - CellMech.a) for ii in eachindex(sol.Density)][idx_for_plots]#[2:4000:end]
t_disc = sol.t[idx_for_plots]

## -----------------------------------------------------------------------------
# Series solution for the first term in dL/dt (continuum limit)
## -----------------------------------------------------------------------------
function L_series(t::Real; P::Int = 200)
    A_inf = (2 / (N * eta_star)) * k_star * ell0_minus_astar
    s = 0.0
    @inbounds for p in 1:P
        s += exp(-α * (2p - 1)^2 * pi^2 * t / (4 * N^2))
    end
    return A_inf * s
end

## -----------------------------------------------------------------------------
# Plotting the series solution and its leading-order short- and long-time asymptotics
## -----------------------------------------------------------------------------

ax2 = Axis(f[1, 2];
    aspect=1,
    xlabel = L"t",
    ylabel = L"|\mathcal{C}(t)| ",
    xscale = log10,
    yscale = log10,
    xminorticksvisible = true,
    yminorticksvisible = true
)

ts_dense = LinRange(1e-3, T, 400)
vals_series = [L_series(t; P=300) for t in ts_dense]
tau = (k, η, N, t) -> (k*pi^2*t) / (4 * N^2 * η)
vals_leading_short = [(((2 / (N * eta_star)) * k_star * ell0_minus_astar) / 4 ) * sqrt(pi / tau(k_star, eta_star, N, t)) for t in ts_dense]
A_inf = (2 / (N * eta_star)) * k_star * ell0_minus_astar
vals_leading_long  = [A_inf * exp(-α * pi^2 * t / (4 * N^2)) for t in ts_dense]

lines!(ax2, ts_dense, max.(abs.(vals_series), 1e-16);
    color = BLACK, linewidth = 5, label = L"\text{Exact form}")
lines!(ax2, ts_dense, abs.(vals_leading_short);
    color = BLUISH_GREEN, linewidth = 5, linestyle = :dash,
    label = L"t \ll \hat{t} \text{ leading order}")
lines!(ax2, ts_dense, abs.(vals_leading_long);
    color = REDDISH_PURPLE, linewidth = 5, linestyle = :dash,
    label = L"t \gg \hat{t} \text{ leading order}")

h(τ) = (sqrt(π)/sqrt(τ)) * erfc(sqrt(τ)/2) - exp(-τ/4)
τ_sol = find_zero(h, 1.0)
# converting τ to t
t_crossover = (2 * N^2 * eta_star * τ_sol) / (k_star * pi^2)
scatter!(ax2, t_crossover, abs.(L_series(t_crossover; P=300));
    color = :grey, marker = :star5, markersize=30, label=L"\hat{t}", strokecolor=BLACK, strokewidth=0.8)
scatter!(ax2, t_disc, abs.(first_term_disc);
    color = YELLOW, markersize=20, label=L"m = %$m", strokecolor=BLACK, strokewidth=0.8)
axislegend(ax2; position = :lb, framevisible=false, patchsize=(40, 10))#, backgroundcolor=:white, framecolor=:white)
display(f)
save(joinpath(FIGDIR, "fig_for_paper_v1.pdf"), f)
