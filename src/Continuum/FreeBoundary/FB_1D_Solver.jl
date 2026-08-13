"""
    FBParams

Parameters for the free boundary continuum model solved in the density `ρ`, on
the fixed computational domain `z ∈ [0,1]` with `x = L(t) z`.

### Fields
- `α`: `k/η`.
- `η`: drag coefficient.
- `k`: spring stiffness.
- `a`: spring resting length.
- `N`: number of spatial grid points.
- `L0`: initial domain length.
- `Δx`: grid spacing on the computational domain, `1/(N-1)`.
- `rBC`, `lBC`: right/left boundary condition, `:fixed` or `:free`.
- `F`, `D`: force and diffusivity as functions of density.
- `P`, `A`: proliferation and apoptosis as functions of density.
- `m`: number of springs per cell.
"""
@kwdef struct FBParams
    α::Float64 = 1.0
    η::Float64 = 1.0
    k::Float64 = 1.0
    a::Float64 = 1.0
    N::Int = 101
    L0::Float64 = 10.0
    Δx::Float64 = 1 / (N - 1)
    rBC::Symbol = :fixed  # right boundary condition (:fixed or :free)
    lBC::Symbol = :fixed  # left boundary condition (:fixed or :free)
    F::Function = ρ -> k * (1 / ρ - a) # Force function from Hookean spring
    D::Function = ρ -> α / ρ^2 # Diffusivity function from Hookean spring
    P::Function = ρ -> 0.0 # Proliferation function
    A::Function = ρ -> 0.0 # Apoptosis function
    m::Int = 100
end

"""
    FBParams_w_correction

As [`FBParams`](@ref), but additionally carrying the initial cell count `Ncells`
and initial density `q₀` needed to evaluate the discrete correction term
[`L_series`](@ref) at the free boundary.

Constructed positionally in field order:
`(α, k, a, η, Δx, L0, N, Ncells, q₀, rBC, F, D, P, A)`.
"""
struct FBParams_w_correction{TF, TD, TP, TA}
    α::Float64
    k::Float64
    a::Float64
    η::Float64
    Δx::Float64
    L0::Float64
    N::Int64
    Ncells::Int64
    q₀::Float64
    rBC::Symbol
    F::TF
    D::TD
    P::TP
    A::TA
end

"""
    rhs_Baker_in_ρ!(du, u, p::FBParams, t)

Spatially discretised right hand side of the continuum free boundary model of
Baker et al. (2019), i.e. without the discrete correction term at the free
boundary.
"""
function rhs_Baker_in_ρ!(du, u, p::FBParams, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    right_BC = p.rBC
    left_BC = p.lBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    m = p.m

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    q_right_ghost = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0
    elseif right_BC == :free
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( (1/η)*F(q[N]))
        dLdt = -(D(q[N])/(q[N] * L))*((q_right_ghost - q[N-1]) / (2*Δx))

    end

    dqidt = 0.0
    for ii in 1:N
        if ii == 1
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (1 / L)*dLdt*( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii]*( P(1/q[ii]) - A(1/q[ii]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

end

"""
    compute_mF_q_term(t, k, η, N, l0, a; t_floor=1e-10, τ=1.1303809932021487)

Leading-order asymptotic approximation to `m F̃(q)` at the free boundary,
switching from the short-time `t^{-1/2}` behaviour to the long-time single-mode
exponential decay at the crossover time `t_switch`.

A small `t_floor` is added to `t` to avoid the `1/√t` singularity at `t = 0`.
"""
@inline function compute_mF_q_term(t::Real, k::Real, η::Real, N::Integer,
                                   l0::Real, a::Real;
                                   t_floor::Real = 1e-10,
                                   τ::Real = 1.1303809932021487)
    # Add the floor to avoid the 1/√t singularity.
    t_eff = t + t_floor

    t_switch = (2 * N^2 * η) / (k * π^2) * τ

    A_tilde = (2/N) * (l0 - a)
    τ_at_t = (k*π^2*t_eff)/(4 * N^2 * η)

    if t_eff < t_switch
        return (A_tilde / 2) * √(pi/τ_at_t)
    else
        return A_tilde * exp(-(k * 0.25 * pi^2 * t_eff)/(2 * N^2 * η))
    end
end

"""
    L_series(t, k, η, N, l0, a; P=200)

Exact eigenfunction series for the discrete correction term `|C(t)|` appearing in
`dL/dt`, truncated at `P` modes:

`(2k/(Nη))(l0 - a) Σ_{p=1}^{P} exp(-α (2p-1)² π² t / (4N²))`,  `α = k/η`.
"""
function L_series(t::Real, k::Real, η::Real, N::Integer, l0::Real, a::Real; P::Int = 200)
    A_inf = (2 / (N * η)) * k * (l0 - a)
    α = k / η
    s = 0.0
    @inbounds for p in 1:P
        s += exp(-α * (2p - 1)^2 * pi^2 * t / (4 * N^2))
    end
    return A_inf * s
end

"""
    rhs_with_correction!(du, u, p::FBParams_w_correction, t)

Spatially discretised right hand side of the continuum free boundary model
*including* the discrete correction term at the free boundary, supplied by
[`L_series`](@ref).

In the `m → ∞` limit the boundary node is a Dirichlet node with `q_N` pinned to
`1/a`, so `du[N]` is frozen at zero rather than integrated.
"""
function rhs_with_correction!(du, u, p::FBParams_w_correction, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    N_IC = p.Ncells
    q₀    = p.q₀
    right_BC = p.rBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    # Initialise before conditional so Julia can infer concrete types
    q_right_ghost = 0.0
    dLdt          = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt          = 0.0

    elseif right_BC == :free
        mF_q_term = L_series(t, k, η, N_IC, 1/q₀, a)
        dqdξ_N    = (q[N] - q[N-1]) / Δx          # one-sided, on the ξ-grid
        dLdt      = -mF_q_term - (D(q[N]) / (2 * q[N] * L)) * dqdξ_N
    end

    dqidt = 0.0

    for ii in 1:N
        if ii == 1
            q_left_ghost = q[2]
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) )
            du[ii] = dqidt

        elseif ii == N
            # m → ∞: Dirichlet boundary, q[N] is pinned to 1/a.
            # Do NOT integrate this node — freeze it so the constraint holds.
            if right_BC == :free
                du[ii] = 0.0
            else
                # :fixed branch keeps the original interior-style update
                if diffusivity_method == "arithmetic"
                    Dm = D(q[ii-1])
                    Di = D(q[ii])
                    Dp = D(q_right_ghost)
                    Dhp = 0.5 * (Di + Dp)
                    Dhm = 0.5 * (Di + Dm)
                elseif diffusivity_method == "harmonic"
                    Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                    Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
                else
                    error("Diffusivity method not recognised")
                end
                z_N = (N - 1) * Δx
                dqidt = (z_N / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                        (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
                du[ii] = 0#dqidt
            end

        else
            z_i = (ii-1) * Δx
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1])
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (mesh-motion / advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                    (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
            du[ii] = dqidt
        end
    end

    du[N + 1] = dLdt

    return nothing
end

"""
    make_initial_condition_FB(N; U0fun = z -> 1.0, L0 = 5.0)

Build the initial state vector `[q(z_1), …, q(z_N), L0]` by sampling `U0fun` on
`N` equally spaced points of the computational domain `z ∈ [0,1]`.
"""
function make_initial_condition_FB(N; U0fun = z -> 1.0 , L0 = 5.0)
    z = range(0.0, 1.0, length=N)
    U0 = [U0fun(zi) for zi in z]
    return vcat(U0, L0)
end
