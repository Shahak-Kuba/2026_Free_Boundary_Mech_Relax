# User defined force function

"""
    F(l, kₛ, a, restoring_force)

Spring restoring force for springs of length `l`, stiffness `kₛ` and resting
length `a`.

- `"hookean"`: `kₛ (l - a)`
- `"nonlinear"`: `kₛ (1/a - 1/l)`
- anything else: `kₛ a² (1/a - 1/l)`
"""
function F(l, kₛ, a, restoring_force)
    if restoring_force == "hookean"
        return (kₛ .* (l' .- a))'
    elseif restoring_force == "nonlinear"
        return (kₛ .* (1 ./ a .- 1 ./ l'))'
    else
        return (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l'))'
    end
end
