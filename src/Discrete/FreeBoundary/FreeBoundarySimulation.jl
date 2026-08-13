"""
    FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)

Run the discrete free boundary simulation from the initial density profile
`FB_IC`, integrating with forward Euler at fixed step `SimTime.δt`.

Returns a [`SimResults_t`](@ref).
"""
function FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions = (false, false))
    # saving callback
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m);
    prob, p = SetupFBODEproblem(FB_IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

    Set_Random_Seed(Seed)

    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)

    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end


"""
    FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)

As [`FreeBoundarySimulation`](@ref), but starting from an explicit vector of
initial node positions `IC` rather than a density profile.
"""
function FreeBoundarySimulation_given_IC(IC::Vector{Float64}, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions = (false, false))
    # saving callback
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m);
    prob, p = SetupFBODEproblem_given_IC(IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

    Set_Random_Seed(Seed)

    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)

    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end
