function _compute_farin_coordinates(nc::NaturalCoordinates{F}, tri::Triangulation, z, ∇) where {F}
    λ = get_coordinates(nc)
    N₀ = get_indices(nc)
    result = is_scalar(z) ? zero(F) : zeros(F, fdim(z))
    for i in eachindex(λ)
        for j in i:lastindex(λ)
            for k in j:lastindex(λ)
                bezier_height, multiplicity = get_contrib(tri, i, j, k, N₀, ∇, z)
                λ_prod = λ[i] * λ[j] * λ[k]
                if is_scalar(z)
                    result += 6bezier_height * λ_prod / multiplicity
                else
                    @. result += 6bezier_height * λ_prod / multiplicity
                end
            end
        end
    end
    return result
end

is_bezier_point(i, j, k) = i == j == k
is_bezier_edge(i, j, k) = (i == j) || (j == k) || (k == i)
is_bezier_face(i, j, k) = (i ≠ j) && (j ≠ k) && (k ≠ i)
function find_bezier_edge(i, j, k)
    if i == j
        return (i, k)
    elseif i == k
        return (i, j)
    else # j == k
        return (j, i)
    end
end
function bezier_point_contribution(i, N₀, z)
    bezier_height = get_data(z, N₀[i])
    return bezier_height, 6
end
function bezier_edge_contribution(tri, i, j, N₀, ∇, z)
    u = N₀[i]
    zᵤ = get_data(z, u)
    bezier_height = zᵤ .+ directional_derivative(tri, i, j, N₀, ∇) ./ 3
    return bezier_height, 2
end
function bezier_face_contribution(tri, i, j, k, N₀, ∇, z)
    u = N₀[i]
    v = N₀[j]
    w = N₀[k]
    point_contrib = (get_data(z, u) .+ get_data(z, v) .+ get_data(z, w)) ./ 3
    edge_contrib = (
        directional_derivative(tri, i, j, N₀, ∇) .+
        directional_derivative(tri, i, k, N₀, ∇) .+
        directional_derivative(tri, j, i, N₀, ∇) .+
        directional_derivative(tri, j, k, N₀, ∇) .+
        directional_derivative(tri, k, i, N₀, ∇) .+
        directional_derivative(tri, k, j, N₀, ∇)
    ) ./ 12
    return point_contrib .+ edge_contrib, 1
end
function get_contrib(tri, i, j, k, N₀, ∇, z)
    if is_bezier_point(i, j, k)
        return bezier_point_contribution(i, N₀, z)
    elseif is_bezier_edge(i, j, k)
        i′, j′ = find_bezier_edge(i, j, k)
        return bezier_edge_contribution(tri, i′, j′, N₀, ∇, z)
    else # is_bezier_face(i,j,k)
        return bezier_face_contribution(tri, i, j, k, N₀, ∇, z)
    end
end

function compute_natural_coordinates(::Farin, tri, interpolation_point, cache=NaturalNeighboursCache(tri); kwargs...)
    return _compute_sibson_coordinates(tri, interpolation_point, cache; kwargs...)
end