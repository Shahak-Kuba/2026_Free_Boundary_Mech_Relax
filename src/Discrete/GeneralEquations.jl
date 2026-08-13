"""
    δ(rᵢ₊₁, rᵢ)

Calculate the Euclidean distance between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The Euclidean distance between the two points.
"""
function δ(rᵢ₊₁, rᵢ)
    val = zeros(1,size(rᵢ,2))
    val .= .√(sum((rᵢ₊₁- rᵢ).^2,dims=1))
    return val
end

"""
    ρ(rᵢ₊₁, rᵢ)

Calculate the reciprocal of the distance (interpreted as density) between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The reciprocal of the distance between the two points.
"""
ρ(rᵢ₊₁, rᵢ) = 1 ./ δ(rᵢ₊₁, rᵢ);

"""
    τ(rᵢ₊₁, rᵢ₋₁)

Calculate the unit tangent vector between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.

# Returns
The unit tangent vector between the two points.
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""
    calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir)

Calculate the spring lengths `l`, densities `ρ`, direction vectors `dv`, tangent
vectors `τ` and normal vectors `n` for the current configuration of interface
nodes.
"""
function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir)
    l = .√(sum((rᵢ₊₁ - rᵢ).^2, dims=1)) # calculating length
    ρ = 1.0 ./ l # calculating density
    dv = (rᵢ₊₁ - rᵢ) ./ l # directional vector
    τ = (rᵢ₊₁ - rᵢ₋₁) ./ .√(sum((rᵢ₊₁ - rᵢ₋₁).^2, dims=1)) # tangent vector

    if growth_dir == "inward" # normal vector
        n = [-1.0, 1.0] .* circshift(dv, 1)
    else
        n = [1.0, -1.0] .* circshift(dv, 1)
    end

    return l, ρ, dv, τ, n # returning all calculated values
end

function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir::ElasticVector{String})
    l = .√(sum((rᵢ₊₁ - rᵢ).^2, dims=1)) # calculating length
    ρ = 1.0 ./ l # calculating density
    dv = (rᵢ₊₁ - rᵢ) ./ l # directional vector
    τ = (rᵢ₊₁ - rᵢ₋₁) ./ .√(sum((rᵢ₊₁ - rᵢ₋₁).^2, dims=1)) # tangent vector

    n = zeros(size(dv))
    for i in eachindex(growth_dir)
        shift_dv = circshift(dv[:, i], 1)
        if growth_dir[i] == "inward"
            n[:, i] = [-1.0, 1.0] .* shift_dv
        else
            n[:, i] = [1.0, -1.0] .* shift_dv
        end
    end

    return l, ρ, dv, τ, n # returning all calculated values
end

"""
    Ω(p)

Calculate the area of a polygon defined by points in `p`.

This function computes the area using the shoelace formula. The polygon is defined by a set of points `p`, and the function iterates through these points to calculate the area.

# Arguments
- `p`: A matrix where each row represents a point of the polygon in 2D space.

# Returns
The absolute area of the polygon.
"""
function Ω(p)
    A = 0
    for ii in axes(p,2)
        if ii == size(p,2)
            A += (p[1,ii]*p[2,1] -  p[2,ii]*p[1,1])
        else
            A += (p[1,ii]*p[2,ii+1] -  p[2,ii]*p[1,ii+1])
        end
    end
    return abs(A)/2;
end
