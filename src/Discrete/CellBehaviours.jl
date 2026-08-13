function Set_Random_Seed(seednum=123)
    Random.seed!(seednum)
end

function calc_cell_lengths(u,m)
    spring_lengths = .√(sum((circshift(u,(0,-1))- u).^2,dims=1))
    return [sum(spring_lengths[i:i+m-1]) for i in 1:m:length(spring_lengths)-m+1]
end

"""
    calc_cell_densities(u, m)

Calculate the cell densities over a window of size `m` for a given state vector `u`.

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `m`: The size of the window over which to calculate the densities.

# Returns
A vector of cell densities, each computed over a window of size `m`.
"""
function calc_cell_densities(u,m)
    return 1 ./ calc_cell_lengths(u,m)
end

"""
    P(event_f, ρ, α, fncs_type)

Calculate the probability of proliferation occurring, given the density `ρ`, the parameter `α`, and the function type `fncs_type`.

If the proliferation event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `α`.
- For `"Constant2"`: `α` only for cells longer than a threshold length.
- For `"Length"`: the probability is proportional to cell length `1/ρ`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the proliferation event is considered to occur.
- `ρ::AbstractVector`: A vector representing densities.
- `α::Float64`: Parameter for scaling the proliferation probability.
- `fncs_type::String`: Type of function used to calculate the probability. Options include:
  - `"Constant"`
  - `"Constant2"`
  - `"Length"`

### Returns
A vector representing the calculated probabilities for proliferation.
"""
function P(event_f, ρ, α, fncs_type)
    if event_f
        if fncs_type == "Constant"
                return α.*ones(size(ρ))

        elseif fncs_type == "Constant2"
                lmin = 10.1
                P_rate = zeros(size(ρ))
                for cell in eachindex(ρ)
                    if 1 ./ ρ[cell] >= lmin
                        P_rate[cell] = α
                    end
                end
                return P_rate

        elseif fncs_type == "Length"
                return α.* (1 ./ ρ)
        end
    else
        return zeros(size(ρ))
    end
end


"""
    A(event_f, ρ, β, fncs_type)

Calculate the probability of apoptosis (cell death) occurring, given the density `ρ`, the parameter `β`, and the function type `fncs_type`.

If the apoptosis event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `β`.
- For `"Constant2"`: `β` only for cells shorter than a threshold length.
- For `"Length"`: the probability increases as the cell shortens below `ld`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the apoptosis event is considered to occur.
- `ρ::AbstractVector`: A vector representing densities.
- `β::Float64`: Parameter for scaling the apoptosis probability.
- `fncs_type::String`: Type of function used to calculate the probability. Options include:
  - `"Constant"`
  - `"Constant2"`
  - `"Length"`

### Returns
A vector representing the calculated probabilities for apoptosis.
"""
function A(event_f, ρ, β, fncs_type)
    if event_f
       if fncs_type == "Constant"
            return β.*ones(size(ρ))
        elseif fncs_type == "Constant2"
            lmax = 10.1
            A_rate = zeros(size(ρ))
            for cell in eachindex(ρ)
                if 1 ./ ρ[cell] <= lmax
                    A_rate[cell] = β
                end
            end
            return A_rate
        elseif fncs_type == "Length"
            ld = 10.1
            A_rate = zeros(size(ρ))
            for ii in eachindex(ρ)
                if 1 ./ ρ[ii] <= ld
                    A_rate[ii] = β*(ld - (1 ./ ρ[ii]))
                else
                    A_rate[ii] = 0
                end
            end
            return A_rate
        end
    else
        return zeros(size(ρ))
    end
end


"""
    E(event_f, ρ, γ, fncs_type)

Calculate the probability of embedding occurring, given the density `ρ`, the parameter `γ`, and the function type `fncs_type`.

If the embedding event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `γ`.
- For `"Constant2"`: `γ` only for cells whose length lies in `[lmin, lmax]`.
- For `"Length"`: `γ` only for cells shorter than `le_max`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the embedding event is considered to occur.
- `ρ::AbstractVector`: A vector representing cell densities.
- `γ::Float64`: Parameter for scaling the embedding probability.
- `fncs_type::String`: Type of function used to calculate the probability as functions of cell length l = 1/q. Options include:
  - `"Constant"`
  - `"Constant2"`
  - `"Length"`

### Returns
A vector representing the calculated probabilities for embedding.
"""
function E(event_f, ρ, γ, fncs_type)
    if event_f
        if fncs_type == "Constant"
                return γ.*ones(size(ρ))
        elseif fncs_type == "Constant2"
                lmin = 0.0
                lmax = 16.1
                E_rate = zeros(size(ρ))
                for cell in eachindex(ρ)
                    if 1 ./ ρ[cell] <= lmax && 1 ./ ρ[cell] >= lmin
                        E_rate[cell] = γ
                    end
                end
                return E_rate
        elseif fncs_type == "Length"
            le_max = 20.0
            le_min = 10.0
            E_rate = zeros(size(ρ))
            for ii in eachindex(ρ)
                if 1 ./ ρ[ii] <= le_max
                    #E_rate[ii] = γ*((le_max - le_min) - (1 ./ ρ[ii]))
                    E_rate[ii] = γ
                else
                    E_rate[ii] = 0.0
                end
            end
        end
    else
        return zeros(size(ρ))
    end
end

function PE(event_f, ρ, γ, fncs_type)
    if event_f
        if fncs_type == "Constant"
            return γ.*ones(size(ρ))
        elseif fncs_type == "Constant2"
            lmin = 0.0
            lmax = 15.1
            PE_rate = zeros(size(ρ))
            for cell in eachindex(ρ)
                if 1 ./ ρ[cell] <= lmax && 1 ./ ρ[cell] >= lmin
                    PE_rate[cell] = γ
                end
            end
            return PE_rate
        end
    else
        return zeros(size(ρ))
    end
end

"""
    cell_probs(uᵢ, curr_t, m, SimTime, prolif, death, embed, prolifembed)

Calculate the probabilities for cell-related events (proliferation, death, embedding,
simultaneous proliferation and embedding) based on the current state vector `uᵢ`.

This function uses `calc_cell_densities` to calculate cell densities and then computes the probabilities for each event. The probabilities are scaled by the time step size `δt`.

### Arguments
- `uᵢ::AbstractVector`: A state vector representing the positions of particles or cells.
- `curr_t::Float64`: The current simulation time.
- `m::Int64`: The size of the window over which to calculate the densities.
- `SimTime::SimTime_t`: An object containing simulation time parameters including event triggering and time step size.
- `prolif::CellEvent_t`: An object representing proliferation event parameters (flag, rate, and function type).
- `death::CellEvent_t`: An object representing death event parameters (flag, rate, and function type).
- `embed::CellEvent_t`: An object representing embedding event parameters (flag, rate, and function type).
- `prolifembed::CellEvent_t`: An object representing simultaneous proliferation/embedding parameters.

### Returns
A tuple of four vectors of probabilities, for proliferation, death, embedding and
simultaneous proliferation/embedding respectively.

If the respective event is not triggered based on the `event_trigger` setting, the function will return vectors of zeros for that event.
"""
function cell_probs(uᵢ,curr_t,m,SimTime,prolif,death,embed,prolifembed)
    ρ = calc_cell_densities(uᵢ,m)
    if SimTime.event_trigger == "Constant"
        return (P(prolif.flag,ρ,prolif.rate,prolif.event_func).*SimTime.event_δt,
        A(death.flag,ρ,death.rate,death.event_func).*SimTime.event_δt,
        E(embed.flag,ρ,embed.rate,embed.event_func).*SimTime.event_δt,
        PE(prolifembed.flag,ρ,prolifembed.rate,prolifembed.event_func).*SimTime.event_δt)
    elseif SimTime.event_trigger == "Periodic"
        if curr_t%SimTime.periodic_δt <= SimTime.event_length
            return (P(prolif.flag,ρ,prolif.rate,prolif.event_func).*SimTime.event_δt,
                A(death.flag,ρ,death.rate,death.event_func).*SimTime.event_δt,
                E(embed.flag,ρ,embed.rate,embed.event_func).*SimTime.event_δt,
                PE(prolifembed.flag,ρ,prolifembed.rate,prolifembed.event_func).*SimTime.event_δt)
        else
            p = zeros(length(P(prolif.flag,ρ,prolif.rate,prolif.event_func)))
            return (p, p, p, p)
        end
    end
end


"""
    find_cell_index(arr::Vector{Float64}, threshold::Float64)

Find the index in `arr` where the cumulative sum first exceeds or equals `threshold`.

This function is typically used in stochastic processes to determine an outcome based on a probability distribution.

# Arguments
- `arr`: A vector of probabilities.
- `threshold`: A threshold value used to find the corresponding index in `arr`.

# Returns
The index in `arr` where the cumulative sum first exceeds or equals `threshold`. If the threshold is not met, the length of `arr` is returned.

"""
function find_cell_index(arr::Vector{Float64}, threshold::Float64)
    cum_sum = cumsum(arr)
    index = findfirst(cum_sum .>= threshold)
    if index === nothing || index > length(arr)
        return index = length(arr)-1  # If the cumulative sum never reaches the threshold
    end
    return index
end


"""
    store_embed_cell_pos(pos)

Stores the position of embedded cells into embedded_cells array.

This function inserts the position vector `pos` of a boundary of an embedded cell into the embedded_cells array.

# Arguments
- `pos`: position vector of the cell that is being embedded into the tissue should contain `m+1` values given `m` springs in the simulation

# Returns
`nothing`. The function modifies `embedded_cells` in place.

!!! note
    `embedded_cells` is a global that is not initialised anywhere in this package.
    Cell embedment is unused by the free boundary mechanical relaxation problem
    (all `CellEvent_t` flags default to `false`), so this path is never reached.
"""
function store_embed_cell_pos(pos)
    global embedded_cells
    insert!(embedded_cells,size(embedded_cells,2),pos)
    return nothing
end

function store_embedded_cell(u, idx, m, t)
    global cell_embedment_times
    push!(cell_embedment_times, t)
    for i = idx:idx+m
        if i > size(u,2)
            store_embed_cell_pos(u[:,1].data)
        else
            store_embed_cell_pos(u[:,i].data)
        end
    end

    return nothing
end

"""
    store_CellMech(u, t, integrator)

`SavingCallback` payload: snapshot the (possibly time-varying) spring mechanical
properties so that post-simulation calculations can use the properties that were
in force at each save time.
"""
function store_CellMech(u,t,integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    data = []
    push!(data,CellMech.kₛ)
    push!(data,CellMech.a)
    push!(data,CellMech.kf)
    push!(data,CellMech.growth_dir)
    return data
end
