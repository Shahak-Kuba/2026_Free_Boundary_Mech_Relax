"""
    FB_ODE!(du, u, p, t)

Right hand side of the discrete free boundary problem. Node 1 is pinned
(`du[:,1] = 0`), interior nodes feel the restoring force of the springs on either
side, and the final node — the free boundary — feels only the spring to its left.
"""
function FB_ODE!(du, u, p , t)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p
    uᵢ₊₁ = circshift(u,(0,-1))
    uᵢ₋₁ = circshift(u,(0,1))

    l,ρ,dv,τ,n = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, "inward")

    du[:,1] = [0. ; 0.]
    du[:,2:end-1] = (1 ./ CellMech.η) .* ( -F(l[1:end-2]',CellMech.kₛ[1:end-2],CellMech.a[1:end-2],CellMech.restoring_force).* dv[:,1:end-2]       +    F(l[2:end-1]',CellMech.kₛ[2:end-1],CellMech.a[2:end-1],CellMech.restoring_force).* dv[:,2:end-1] )
    du[:,end] = (1 ./ CellMech.η) .* ( -F(l[end-1],CellMech.kₛ[end-1],CellMech.a[end-1],CellMech.restoring_force).* dv[:,end-1] )

end
