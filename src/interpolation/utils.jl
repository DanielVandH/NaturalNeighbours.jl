function compute_bowyer_envelope!(envelope, tri::Triangulation, history::InsertionEventHistory, temp_adjacent::Adjacent, point; kwargs...) #kwargs are add_point! kwargs
    empty!(history)
    empty!(envelope)
    empty!(get_adjacent(temp_adjacent))
    n = num_points(tri)
    I = integer_type(tri)
    V = add_point!(tri, point; store_event_history=Val(true), event_history=history, peek=Val(true), kwargs...)
    all_triangles = each_added_triangle(history)
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

function polygon_area(points)
    n = num_points(points)
    p, q, r, s = get_point(points, 1, 2, n, n - 1)
    px, py = getxy(p)
    _, qy = getxy(q)
    rx, ry = getxy(r)
    _, sy = getxy(s)
    area = px * (qy - ry) + rx * (py - sy)
    for i in 2:(n-1)
        p, q, r = get_point(points, i, i + 1, i - 1)
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

function handle_duplicate_points!(tri, interpolation_point, coordinates::AbstractVector{F}, envelope) where {F}
    idx = findfirst(i -> get_point(tri, i) == getxy(interpolation_point), envelope)
    envelope_idx = envelope[idx]
    resize!(coordinates, 1)
    resize!(envelope, 1)
    envelope[begin] = envelope_idx
    coordinates[begin] = one(F)
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end