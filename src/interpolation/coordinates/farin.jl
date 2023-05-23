function _compute_farin_coordinates(nc::NaturalCoordinates{F}, tri::Triangulation, z, ∇) where {F}
    λ = get_coordinates(nc)
    N₀ = get_indices(nc)
    result = zero(F)
    for i in eachindex(λ)
        for j in i:lastindex(λ)
            for k in j:lastindex(λ)
                bezier_height, num_occurrences = get_contrib(tri, i, j, k, N₀, ∇, z)
                λ_prod = λ[i] * λ[j] * λ[k]
                scale = 6bezier_height / num_occurrences
                result += scale * λ_prod
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
    bezier_height = z[N₀[i]]
    num_occurrences = 6
    return bezier_height, num_occurrences
end
function bezier_edge_contribution(tri, i, j, N₀, ∇, z)
    u = N₀[i]
    v = N₀[j]
    p, q = get_point(tri, u, v)
    px, py = getxy(p)
    qx, qy = getxy(q)
    dx = qx - px
    dy = qy - py
    ∇ᵤ = ∇[u]
    ∇ᵤx, ∇ᵤy = getxy(∇ᵤ)
    zᵤ = z[u]
    bezier_height = zᵤ + (1 / 3) * (dx * ∇ᵤx + dy * ∇ᵤy)
    num_occurrences = 2
    return bezier_height, num_occurrences
end
function bezier_face_contribution(tri, i, j, k, N₀, ∇, z)
    edge_contrib = 0.0
    for (u, v) in ((i, j), (i, k), (j, i), (j, k), (k, i), (k, j))
        edge_contrib += bezier_edge_contribution(tri, u, v, N₀, ∇, z)[1]
    end
    edge_contrib /= 4
    u = N₀[i]
    v = N₀[j]
    w = N₀[k]
    point_contrib = (1 / 6) * (z[u] + z[v] + z[w])
    bezier_height = edge_contrib - point_contrib
    num_occurrences = 1
    return bezier_height, num_occurrences
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