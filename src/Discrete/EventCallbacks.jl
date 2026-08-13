# Cell behaviour Callback

"""
    event_affect!(integrator)

Update the state of `integrator` based on probabilistic cellular events.

This function modifies `integrator` in place. It uses the parameters and state from `integrator` to compute probabilities for different cellular events: proliferation (`prolif`), death (`death`), and embedding (`embed`). Based on these probabilities and random draws, it updates the state of `integrator`.

# Arguments
- `integrator`: The integrator object containing the current state and parameters.

# Returns
`nothing`. The function modifies `integrator` in place.

!!! note
    For the free boundary mechanical relaxation problem all `CellEvent_t` flags are
    `false`, so `cell_probs` returns zeros and the body below never fires. The
    callback is retained so that the solve path is unchanged from the original
    `MechCellTissueGrowth.jl` implementation.
"""
function event_affect!(integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    u = integrator.u
    curr_t = integrator.t
    (p, a, e, pe) = cell_probs(u, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    (r1, r2, r3) = rand(3)
    total_prob = sum(p) + sum(a) + sum(e) + sum(pe)

    if r1 < total_prob # check if event has occurred
        event_type = r2 * total_prob
        event_occuring = ""
        if event_type < sum(p) # prolif occurs
            # finding proliferating cell index
            idx = find_cell_index(p, r3 * sum(p))
            left_spring_index = idx * Domain.m - (Domain.m - 1)
            right_spring_index = idx * Domain.m + 1
            right_spring_index = right_spring_index > size(u, 2) ? 1 : right_spring_index
            # interpolating & inserting new daughter cell spring boundaries
            centre_spring_pos = Domain.m % 2 == 0 ? u[:, Int64(left_spring_index + Domain.m / 2)] : (u[:, Int64(left_spring_index + (Domain.m - 1) / 2)] + u[:, Int64(left_spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            interp_points_left = LinearInterp(u[:, left_spring_index], centre_spring_pos, Domain.m)
            interp_points_right = LinearInterp(centre_spring_pos, u[:, right_spring_index], Domain.m)
            for i = 1:Domain.m - 1
                deleteat!(u, left_spring_index + 1)
            end
            new_prolif_spring_pos = vcat(interp_points_left, [centre_spring_pos], interp_points_right)
            for point in reverse(new_prolif_spring_pos)
                insert!(u, left_spring_index + 1, point)
            end
            # inserting new daughter cell spring mechanical properties
            for i in 1:Domain.m
                insert!(CellMech.kₛ, left_spring_index + 1, CellMech.kₛ[idx])
                insert!(CellMech.kf, left_spring_index + 1, CellMech.kf[idx])
                insert!(CellMech.a, left_spring_index + 1, CellMech.a[idx])
                insert!(CellMech.growth_dir, left_spring_index + 1, CellMech.growth_dir[idx])
            end
            event_occuring = "prolif"
        elseif event_type < sum(p) + sum(a) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            spring_index = idx * Domain.m - (Domain.m - 1)
            boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            if spring_index > size(u, 2)
                u[:, 1] .= boundaries_midpoint
            else
                u[:, spring_index] .= boundaries_midpoint
            end
            event_occuring = "death"

        elseif event_type < sum(p) + sum(a) + sum(e) # embed only occurs
            idx = find_cell_index(e, r3 * sum(e))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
            boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            if spring_index > size(u, 2)
                u[:, 1] .= boundaries_midpoint
            else
                u[:, spring_index] .= boundaries_midpoint
            end
            event_occuring = "embed"
        else # prolif and embed simultaneously occur
            idx = find_cell_index(pe, r3 * sum(pe))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
        end
        resize!(integrator, (2, size(integrator.u, 2)))

        if size(integrator.u, 2) != size(CellMech.kₛ, 1)
            error("sizing error with: ", event_occuring)
        end
    end
    nothing
end

"""
    LinearInterp(a, b, m)

Return the `m-1` interior points that linearly interpolate between `a` and `b`.
"""
function LinearInterp(a::ElasticVector, b::ElasticVector, m::Int64)
    dt = 1 / m
    T = dt:dt:1-dt
    Interp_P = [a + t * (b - a) for t in T]
    return Interp_P
end

terminate_affect!(integrator) = terminate!(integrator)
