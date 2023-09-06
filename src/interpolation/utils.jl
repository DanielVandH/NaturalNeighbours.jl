function compute_bowyer_envelope!(envelope, tri::Triangulation, history::InsertionEventHistory, temp_adjacent::Adjacent, point; kwargs...) #kwargs are add_point! kwargs
    empty!(history)
    empty!(envelope)
    empty!(get_adjacent(temp_adjacent))
    n = num_points(tri)
    I = integer_type(tri)
    V = add_point!(tri, point; store_event_history=Val(true), event_history=history, peek=Val(true), kwargs...)
    all_triangles = each_added_triangle(history)
    if isempty(all_triangles) # This is possible for constrained triangulations
        return envelope, temp_adjacent, history, V
    end
    for T in all_triangles
        add_triangle!(temp_adjacent, T)
    end
    T = first(all_triangles)
    i, j, _ = indices(T)
    v = i == I(n + 1) ? j : i # get a vertex on the envelope
    push!(envelope, v)
    for i in 2:num_triangles(all_triangles)
        v = get_adjacent(temp_adjacent, I(n + 1), v)
        push!(envelope, v)
    end
    push!(envelope, envelope[begin])
    return envelope, temp_adjacent, history, V
end
function compute_bowyer_envelope!(envelope, tri::Triangulation, point; kwargs...)
    I = integer_type(tri)
    E = edge_type(tri)
    A = Adjacent{I,E}()
    return compute_bowyer_envelope!(envelope, tri, initialise_event_history(tri), A, point; kwargs...)
end
function compute_bowyer_envelope(tri::Triangulation, point; kwargs...)
    I = integer_type(tri)
    envelope = I[]
    return compute_bowyer_envelope!(envelope, tri, point; kwargs...)
end

_num_points(points::Tuple) = length(points)
_num_points(points) = num_points(points)
_get_point(points::Tuple, i...) = ntuple(j->points[i[j]], length(i))
_get_point(points, i...) = get_point(points, i...)
function polygon_area(points)
    n = _num_points(points)
    p, q, r, s = _get_point(points, 1, 2, n, n - 1)
    px, py = getxy(p)
    _, qy = getxy(q)
    rx, ry = getxy(r)
    _, sy = getxy(s)
    area = px * (qy - ry) + rx * (py - sy)
    for i in 2:(n-1)
        p, q, r = _get_point(points, i, i + 1, i - 1)
        px, py = getxy(p)
        _, qy = getxy(q)
        rx, ry = getxy(r)
        area += px * (qy - ry)
    end
    return area / 2
end

function get_barycentric_deviation(natural_coordinates::NaturalCoordinates{F}) where {F}
    coordinates = get_coordinates(natural_coordinates)
    indices = get_indices(natural_coordinates)
    interpolation_point = get_interpolation_point(natural_coordinates)
    triangulation = get_triangulation(natural_coordinates)
    x̂ = zero(F)
    ŷ = zero(F)
    for (λ, k) in zip(coordinates, indices)
        p = get_point(triangulation, k)
        px, py = getxy(p)
        x̂ += λ * px
        ŷ += λ * py
    end
    x, y = getxy(interpolation_point)
    δ² = (x - x̂)^2 + (y - ŷ)^2
    return sqrt(δ²)
end

function handle_duplicate_points!(tri, interpolation_point, coordinates::AbstractVector{F}, envelope, u, prev_u, next_u) where {F}
    p, q, r = get_point(tri, u, prev_u, next_u)
    xy = getxy(interpolation_point)
    envelope_idx = if xy == p
        u
    elseif xy == q
        prev_u
    elseif xy == r
        next_u
    else
        idx = findfirst(i -> get_point(tri, i) == getxy(interpolation_point), envelope)
        envelope[idx]
    end
    resize!(coordinates, 1)
    resize!(envelope, 1)
    envelope[begin] = envelope_idx
    coordinates[begin] = one(F)
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end

function sort_five(i, j, k, ℓ, m)
    if j < i
        i, j = j, i
    end
    if k < i
        i, k = k, i
    end
    if ℓ < i
        i, ℓ = ℓ, i
    end
    if m < i
        i, m = m, i
    end
    if k < j
        j, k = k, j
    end
    if ℓ < j
        j, ℓ = ℓ, j
    end
    if m < j
        j, m = m, j
    end
    if ℓ < k
        k, ℓ = ℓ, k
    end
    if m < k
        k, m = m, k
    end
    if m < ℓ
        ℓ, m = m, ℓ
    end
    return i, j, k, ℓ, m
end

function count_unique_sorted(i, j, k, ℓ, m) # assumes sorted
    n = 5
    if i == j
        n -= 1
    end
    if j == k
        n -= 1
    end
    if k == ℓ
        n -= 1
    end
    if ℓ == m
        n -= 1
    end
    return n
end

#=
Returns (standard_sort, case)
Standard forms:
    Case 1. iiiii 
    Case 2. iiiij 
    Case 3. iiijj 
    Case 4. iiijk 
    Case 5. iijjk 
    Case 6. iijkℓ 
    Case 7. ijkℓm
=#
function group_sort(i, j, k, ℓ, m)
    i′, j′, k′, ℓ′, m′ = sort_five(i, j, k, ℓ, m)
    num_unique = count_unique_sorted(i′, j′, k′, ℓ′, m′)
    if num_unique == 1
        # uuuuu 
        u = i′
        return (u, u, u, u, u), 1
    elseif num_unique == 2
        if i′ == j′ == k′ == ℓ′         # iiiij
            u, v = i′, m′
            return (u, u, u, u, v), 2
        elseif j′ == k′ == ℓ′ == m′     # jiiii
            u, v = j′, i′
            return (u, u, u, u, v), 2
        elseif i′ == j′ == k′           # iiijj
            u, v = i′, ℓ′
            return (u, u, u, v, v), 3
        else                            # jjiii
            u, v = k′, i′
            return (u, u, u, v, v), 3
        end
    elseif num_unique == 3
        if i′ == j′ == k′               # iiijk
            u, v, w = i′, ℓ′, m′
            return (u, u, u, v, w), 4
        elseif k′ == ℓ′ == m′           # jkiii
            u, v, w, = k′, i′, j′
            return (u, u, u, v, w), 4
        elseif j′ == k′ == ℓ′           # ijjjk 
            u, v, w = j′, i′, m′
            return (u, u, u, v, w), 4
        elseif (i′ == j′) && (k′ == ℓ′) # iijjk
            u, v, w = i′, k′, m′
            return (u, u, v, v, w), 5
        elseif (j′ == k′) && (ℓ′ == m′) # kiijj
            u, v, w = j′, ℓ′, i′
            return (u, u, v, v, w), 5
        else                            # iikjj
            u, v, w = i′, ℓ′, k′
            return (u, u, v, v, w), 5
        end
    elseif num_unique == 4
        if i′ == j′                     # iijkℓ
            u, v, w, x = i′, k′, ℓ′, m′
        elseif j′ == k′                 # jiikℓ
            u, v, w, x = j′, i′, ℓ′, m′
        elseif k′ == ℓ′                 # jkiiℓ
            u, v, w, x = k′, i′, j′, m′
        else                            # jkℓii
            u, v, w, x = ℓ′, i′, j′, k′
        end
        return (u, u, v, w, x), 6
    else                                # ijkℓm
        return (i′, j′, k′, ℓ′, m′), 7
    end
end

# computes dot(∇ᵢ, xⱼ - xᵢ)
function directional_derivative(tri, i, j, N₀, ∇) # zᵢⱼ
    u = N₀[i]
    v = N₀[j]
    p, q = get_point(tri, u, v)
    px, py = getxy(p)
    qx, qy = getxy(q)
    dx = qx - px
    dy = qy - py
    ∇ᵤ = ∇[u]
    ∇ᵤx, ∇ᵤy = getxy(∇ᵤ)
    return dx * ∇ᵤx + dy * ∇ᵤy
end

# computes dot(xⱼ - xᵢ, B(xₖ - xᵢ))
function hessian_form(tri, u, v, w, N₀, H) # zᵢ,ⱼₖ
    i = N₀[u]
    j = N₀[v]
    k = N₀[w]
    xᵢ, xⱼ, xₖ = get_point(tri, i, j, k)
    B = H[i]
    B₁₁, B₂₂, B₁₂ = B
    dxᵢⱼ = xⱼ[1] - xᵢ[1]
    dyᵢⱼ = xⱼ[2] - xᵢ[2]
    dxᵢₖ = xₖ[1] - xᵢ[1]
    dyᵢₖ = xₖ[2] - xᵢ[2]
    return dxᵢⱼ * (dxᵢₖ * B₁₁ + dyᵢₖ * B₁₂) + dyᵢⱼ * (dxᵢₖ * B₁₂ + dyᵢₖ * B₂₂)
end