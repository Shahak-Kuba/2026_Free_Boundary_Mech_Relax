
"""
    SimResults_t

A data structure representing the results of a simulation.

This structure contains various fields that store the results and relevant data from a simulation process. It includes the type of simulation, time steps, state vectors, forces, densities, and other relevant quantities.

# Fields
- `btype`: A `String` indicating the type of simulation or boundary condition.
- `t`: A `Vector{Float64}` representing the time steps at which the simulation results are recorded.
- `u`: A `Vector` of `Matrix{Float64}` representing the state of the system at each time step.
- `∑F`: A `Vector` of `Vector{Float64}`, representing the sum of forces at each time step.
- `Density`: A `Vector` of `Vector{Float64}` representing the density of the system at each time step.
- `Vₙ`: A `Vector{Vector{Float64}}` representing the normal velocity at each time step.
- `Ω`: A `Vector{Float64}` representing the void area within the pore at each time step.
- `σ`: A `Vector` of `Vector{Float64}`, representing the stress of the springs/ cells at each time step.
- `Κ`: A `Vector` of `Vector{Float64}`, representing the approximate curvature of the moving boundary at each time step.
- `CellCount`: A `Vector{Int64}` representing the amount of active cells at a give time `t`.

"""
struct SimResults_t
    btype::String
    t::Vector{Float64}
    u::Vector{Matrix{Float64}}
    ∑F::Vector{Vector{Float64}}
    Density::Vector{Vector{Float64}}
    Vₙ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    σ::Vector{Vector{Float64}}
    Κ::Vector{Vector{Float64}}
    CellCount::Vector{Int64}
end

"""
    CellEvent_t

Encapsulates the properties related to cell events in a simulation.
Includes parameters that define whether an event is active, the rate
of the event, and the function associated with the event.

### Fields
- `flag::Bool`: Indicates whether the cell event is active. (Default: `false`)
- `rate::Float64`: The rate at which the event occurs. (Default: `0.0`)
- `event_func::String`: The function that defines the event's behavior. (Default: `"Constant"`)

"""
@kwdef struct CellEvent_t
    flag::Bool = false
    rate::Float64 = 0.0
    event_func::String = "Constant"
   #event_dep::String = "Density"
end

"""
    CellMechProperties_t

Holds the mechanical properties of a cell, including spring constants,
damping coefficients, and growth directions.

### Fields
- `kₛ::Float64`: Spring constant representing the stiffness of the cell. (Default: `25.0`)
- `η::Float64`: Damping coefficient for the mechanical response. (Default: `1.0`)
- `restoring_force::String`: Defines the nature of the restoring force. (Default: `"nonlinear"`)
- `kf::Float64`: A constant related to the restoring force calculation. (Default: `87.842`)
- `growth_dir::String`: The direction of cell growth. (Default: `"inward"`)
- `D::Float64`: Diffusion coefficient or related parameter. (Default: `0.0`)
- `a::Float64`: A parameter that might relate to cell size or other mechanical properties. (Default: `10.0`)

"""
@kwdef struct CellMechProperties_t
    kₛ::Float64 = 25.0
    η::Float64 = 1.0
    restoring_force::String = "nonlinear"
    kf::Float64 = 87.842
    growth_dir::String = "inward"
    D::Float64 = 0.0
    a::Float64 = 10.0
end

@kwdef mutable struct HeterogeneousCellMechProperties_t
    kₛ::ElasticArray{Float64}
    η::Float64 = 1.0
    restoring_force::String = "nonlinear"
    kf::ElasticArray{Float64}
    growth_dir::ElasticArray{String}
    a::ElasticArray{Float64}
end

"""
    DomainProperties_t

Defines the properties of the simulation domain, including dimensions,
shape, and type.

### Fields
- `N::Int64`: Total number of cells in the simulation. (Default: `100`)
- `m::Int64`: Number of springs connected to each cell. (Default: `6`)
- `R₀::Float64`: The shape radius in micrometers (μm). (Default: `105.0`)
- `domain_type::String`: Type of simulation domain. (Default: `"2D"`)
- `btype::String`: Shape of the domain. (Default: `"circle"`)
  - **Options:** `["circle", "triangle", "square", "hex", "star", "cross"]`
- `dist_type::String`: Distribution type for cells within the domain. (Default: `"Linear"`)

"""
@kwdef struct DomainProperties_t
    N::Int64 = 100
    m::Int64 = 6 # number of springs per cell
    R₀::Float64 = 105.0  # shape radius μm
    domain_type::String = "2D"
    btype::String = "circle"  #Options: ["circle", "triangle", "square", "hex", "star","cross"]
    dist_type::String = "Linear"
    Ω_0::Float64 = π*R₀^2
    dir_to_EXP_Image::String = ""
end

"""
    SimTime_t

Specifies the time-related properties for the simulation, including
total time, time step, and event trigger settings.

### Fields
- `Tmax::Float64`: Maximum simulation time in days. (Default: `20.0`)
- `δt::Float64`: Time step for the simulation. (Default: `0.01`)
- `event_trigger::String`: The method for triggering events during the simulation. (Default: `"Constant"`)
- `event_δt::Float64`: Time interval for triggering events. (Default: `0.01`)
- `event_length::Float64`: Duration of the event. Should be greater than 0 if `event_trigger` is set to "Periodic". (Default: `0.0`)

"""
@kwdef struct SimTime_t
    Tmax::Float64 = 20.0 # days
    δt::Float64 = 0.01
    event_trigger::String = "Constant"
    event_δt::Float64 = 0.01
    periodic_δt::Float64 = 4.0
    event_length:: Float64 = 0.0 # make this > 0 in the case that event_trigger is Periodic
end

"""
    FB_IC_t

Initial condition for the free boundary problem: an initial cell density profile
`q0` on `[0, L0]`, its derivative `q0_der`, and the initial domain length `L0`.
"""
@kwdef struct FB_IC_t
    q0::Function = x -> 2.0
    q0_der::Function = x -> 0.0
    L0::Float64 = 10.0
end
