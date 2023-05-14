module NaturalNeighbourInterp

import DelaunayTriangulation: DelaunayTriangulation,
    triangulate,
    integer_type,
    num_points,
    InsertionEventHistory,
    add_point!,
    each_added_triangle,
    indices,
    initialise_event_history,
    Triangulation,
    triangulate,
    is_boundary_index,
    construct_triangle,
    Adjacent,
    add_triangle!,
    edge_type,
    get_adjacent,
    num_triangles,
    number_type,
    previndex_circular,
    nextindex_circular,
    get_point,
    triangle_circumcenter,
    num_points,
    number_type,
    getxy,
    polygon_features,
    getpoint,
    num_solid_vertices,
    has_ghost_triangles,
    add_ghost_triangles!,
    is_collinear,
    has_boundary_nodes,
    get_triangulation,
    rotate_ghost_triangle_to_standard_form,
    is_on,
    is_degenerate,
    point_position_relative_to_triangle,
    point_position_on_line_segment,
    point_position_relative_to_line,
    is_ghost_triangle,
    distance_to_polygon,
    initial,
    terminal,
    find_edge,
    triangle_type,
    is_boundary_triangle,
    replace_boundary_triangle_with_ghost_triangle,
    each_solid_triangle,
    jump_and_march
import ChunkSplitters: chunks

num_points(::NTuple{N,F}) where {N,F} = N
getpoint(p::NTuple{N,F}, i::Integer) where {N,F} = p[i]

export interpolate

struct NaturalCoordinates{F,I,T<:Triangulation}
    coordinates::Vector{F}
    indices::Vector{I}
    interpolation_point::NTuple{2,F}
    triangulation::T
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalCoordinates{F,I}) where {F,I}
    coordinates = get_coordinates(nc)
    indices = get_indices(nc)
    interpolation_point = get_interpolation_point(nc)
    println(io, "NaturalCoordinates{", F, ",", I, "}")
    println(io, "    u: ", interpolation_point)
    println(io, "    λ: ", coordinates)
    print(io, "    k: ", indices)
end
get_coordinates(nc::NaturalCoordinates) = nc.coordinates
get_indices(nc::NaturalCoordinates) = nc.indices
get_interpolation_point(nc::NaturalCoordinates) = nc.interpolation_point
get_triangulation(nc::NaturalCoordinates) = nc.triangulation

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

# pre-insertion area component from the envelope[i]th generator
function pre_insertion_area!(poly_points, envelope, i, tri::Triangulation)
    empty!(poly_points)
    u = envelope[i]
    prev_u = envelope[previndex_circular(envelope, i)]
    next_u = envelope[nextindex_circular(envelope, i)]
    v = next_u
    ux, uy = get_point(tri, u)
    vx, vy = get_point(tri, v)
    mx1, my1 = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx1, my1))
    while v ≠ prev_u
        w = get_adjacent(tri, u, v)
        cx, cy = triangle_circumcenter(tri, (u, v, w))
        push!(poly_points, (cx, cy))
        v = w
    end
    vx, vy = get_point(tri, v)
    mx, my = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx, my))
    push!(poly_points, (mx1, my1))
    return polygon_area(poly_points)
end

# post-insertion area component from the envelope[i]th generator
function post_insertion_area(envelope, i, tri::Triangulation, interpolation_point)
    u = envelope[i]
    prev_u = envelope[previndex_circular(envelope, i)]
    next_u = envelope[nextindex_circular(envelope, i)]
    p, q, r = get_point(tri, u, prev_u, next_u)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    mpq = (px + qx) / 2, (py + qy) / 2
    mpr = (px + rx) / 2, (py + ry) / 2
    g1 = triangle_circumcenter(p, r, interpolation_point)
    F = number_type(tri)
    if any(isnan, g1)
        # The circumcenter is NaN when the triangle is degenerate, 
        # meaning one of the points is a duplicate with another.
        # Since the triangulation is assumed to be valid, it must be that
        # interpolation_point is one of p or r. In particular, the new point 
        # is just one of the others, and so there will be no changes 
        # in the area. We return NaN as a flag.
        return F(NaN)
    end
    g2 = triangle_circumcenter(q, p, interpolation_point)
    any(isnan, g2) && return F(NaN)
    points = (mpq, mpr, g1, g2, mpq)
    return polygon_area(points)
end

struct InterpolantCache{F,I,H,E,R}
    coordinates::Vector{F}
    envelope::Vector{I}
    insertion_event_history::H
    poly_points::Vector{NTuple{2,F}}
    temp_adjacent::Adjacent{I,E}
    last_triangle::R
end
get_coordinates(cache::InterpolantCache) = cache.coordinates
get_envelope(cache::InterpolantCache) = cache.envelope
get_insertion_event_history(cache::InterpolantCache) = cache.insertion_event_history
get_poly_points(cache::InterpolantCache) = cache.poly_points
get_temp_adjacent(cache::InterpolantCache) = cache.temp_adjacent
get_last_triangle(cache::InterpolantCache) = cache.last_triangle
function InterpolantCache(tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}
    coordinates = number_type(tri)[]
    envelope = I[]
    insertion_event_history = initialise_event_history(tri)
    poly_points = NTuple{2,number_type(tri)}[]
    temp_adjacent = Adjacent{I,E}()
    last_triangle = (Ref ∘ indices ∘ first ∘ each_solid_triangle)(tri)
    return InterpolantCache(coordinates, envelope, insertion_event_history, poly_points, temp_adjacent, last_triangle)
end

function two_point_interpolate!(tri, i, j, c)
    # Project c onto the line through (a, b): https://stackoverflow.com/a/15187473
    # The orthogonal projection is not guaranteed to be on the line segment (corresponding 
    # to t < 0 or t > 1), in which case the weights are no longer a convex combination.
    a, b = get_point(tri, i, j)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    ℓ² = (ax - bx)^2 + (ay - by)^2
    t = (cx - ax) * (bx - ax) + (cy - ay) * (by - ay)
    t /= ℓ²
    return t
end

function two_point_interpolate!(coordinates, envelope, tri, i, j, r) #interpolate r using two points i, j
    t = two_point_interpolate!(tri, i, j, r)
    resize!(coordinates, 2)
    resize!(envelope, 2)
    coordinates[1] = one(t) - t
    coordinates[2] = t
    envelope[1] = i
    envelope[2] = j
    return NaturalCoordinates(coordinates, envelope, r, tri)
end

function compute_natural_coordinates(
    tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL},
    interpolation_point,
    cache::InterpolantCache{F}=InterpolantCache(tri);
    method=:sibson,
    kwargs...
) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,F}
    if method == :sibson
        return _compute_sibson_coordinates(tri, interpolation_point, cache; kwargs...)
    elseif method == :triangle
        return _compute_triangle_coordinates(tri, interpolation_point, cache; kwargs...)
    else
        throw(ArgumentError("method must be one of :sibson or :triangle."))
    end
end

function check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    if is_ghost_triangle(V)
        V = rotate_ghost_triangle_to_standard_form(V)
        i, j, _ = indices(V)
        last_triangle[] = (i, j, get_adjacent(tri, j, i))
    else
        last_triangle[] = indices(V)
    end
    if is_boundary_triangle(tri, V)
        _V = replace_boundary_triangle_with_ghost_triangle(tri, V)
        _u, _w, _ = indices(_V)
        cert = point_position_relative_to_line(tri, _u, _w, interpolation_point)
        if is_collinear(cert)
            cert = point_position_on_line_segment(tri, _u, _w, interpolation_point)
            if is_on(cert) || is_degenerate(cert)
                V = _V
            end
        end
    end
    if is_ghost_triangle(V)
        V = rotate_ghost_triangle_to_standard_form(V)
        i, j, _ = indices(V)
        return i, j, true
    end
    i, j, _ = indices(V)
    return i, j, false
end

function _compute_sibson_coordinates(
    tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL},
    interpolation_point,
    cache::InterpolantCache{F}=InterpolantCache(tri);
    kwargs...
) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,F}
    coordinates = get_coordinates(cache)
    envelope = get_envelope(cache)
    insertion_event_history = get_insertion_event_history(cache)
    poly_points = get_poly_points(cache)
    temp_adjacent = get_temp_adjacent(cache)
    last_triangle = get_last_triangle(cache)
    envelope, temp_adjacent, insertion_event_history, V = compute_bowyer_envelope!(envelope, tri, insertion_event_history, temp_adjacent, interpolation_point; try_points=last_triangle[], kwargs...) #kwargs are add_point! kwargs
    i, j, return_flag = check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    return_flag && return two_point_interpolate!(coordinates, envelope, tri, i, j, interpolation_point)
    resize!(coordinates, length(envelope) - 1)
    w = zero(number_type(tri))
    duplicate_point_flag = false
    for i in firstindex(envelope):(lastindex(envelope)-1)
        pre = pre_insertion_area!(poly_points, envelope, i, tri)
        post = post_insertion_area(envelope, i, tri, interpolation_point)
        if isnan(post)
            duplicate_point_flag = true
            break
        end
        coordinates[i] = pre - post
        w += coordinates[i]
    end
    pop!(envelope)
    if duplicate_point_flag
        idx = findfirst(i -> get_point(tri, i) == getxy(interpolation_point), envelope)
        fill!(coordinates, zero(F))
        coordinates[idx] = one(F)
    else
        coordinates ./= w
    end
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end

function _compute_triangle_coordinates(
    tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL},
    interpolation_point,
    cache::InterpolantCache{F}=InterpolantCache(tri);
    kwargs...
) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,F}
    coordinates = get_coordinates(cache)
    envelope = get_envelope(cache)
    last_triangle = get_last_triangle(cache)
    V = jump_and_march(tri, interpolation_point; try_points=last_triangle[])
    i, j, return_flag = check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    return_flag && return two_point_interpolate!(coordinates, envelope, tri, i, j, interpolation_point)
    i, j, k = indices(V)
    resize!(coordinates, 3)
    resize!(envelope, 3)
    p, q, r = get_point(tri, i, j, k)
    x₁, y₁ = getxy(p)
    x₂, y₂ = getxy(q)
    x₃, y₃ = getxy(r)
    x, y = getxy(interpolation_point)
    Δ = (y₂ - y₃) * (x₁ - x₃) + (x₃ - x₂) * (y₁ - y₃)
    λ₁ = ((y₂ - y₃) * (x - x₃) + (x₃ - x₂) * (y - y₃)) / Δ
    λ₂ = ((y₃ - y₁) * (x - x₃) + (x₁ - x₃) * (y - y₃)) / Δ
    λ₃ = one(λ₁) - λ₁ - λ₂
    coordinates[1] = λ₁
    coordinates[2] = λ₂
    coordinates[3] = λ₃
    envelope[1] = i
    envelope[2] = j
    envelope[3] = k
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
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

struct NaturalNeighbourInterpolant{T<:Triangulation,F,C}
    triangulation::T
    z::Vector{F}
    cache::C
    function NaturalNeighbourInterpolant(tri::T, z::Vector{F}) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        caches = [InterpolantCache(tri) for _ in 1:nt]
        return new{T,F,typeof(caches)}(tri, z, caches)
    end
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalNeighbourInterpolant)
    z = get_z(nc)
    println(io, "Natural Neighbour Interpolant")
    print(io, "    z: ", z)
end
get_triangulation(ni::NaturalNeighbourInterpolant) = ni.triangulation
get_z(ni::NaturalNeighbourInterpolant) = ni.z
get_cache(ni::NaturalNeighbourInterpolant) = ni.cache
get_cache(ni::NaturalNeighbourInterpolant, id) = ni.cache[id]
interpolate(tri::Triangulation, z) = NaturalNeighbourInterpolant(tri, z)
function interpolate(points, z)
    tri = triangulate(points, delete_ghosts=false)
    return interpolate(tri, z)
end
function interpolate(x::AbstractVector, y::AbstractVector, z)
    @assert length(x) == length(y) == length(z) "x, y, and z must have the same length."
    points = [(ξ, η) for (ξ, η) in zip(x, y)]
    return interpolate(points, z)
end

function _eval_interp(itp::NaturalNeighbourInterpolant, p, cache; method=:sibson, kwargs...)
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(tri, p, cache; method, kwargs...)
    z = get_z(itp)
    coordinates = get_coordinates(nc)
    indices = get_indices(nc)
    F = number_type(tri)
    val = zero(F)
    for (λ, k) in zip(coordinates, indices)
        zₖ = z[k]
        val += λ * zₖ
    end
    return val
end

function (itp::NaturalNeighbourInterpolant)(x, y, id::Integer=1; parallel=false, method=:sibson, kwargs...)
    p = (x, y)
    cache = get_cache(itp, id)
    return _eval_interp(itp, p, cache; method, kwargs...)
end

function (itp::NaturalNeighbourInterpolant)(vals::AbstractVector, x::AbstractVector, y::AbstractVector; parallel=true, method=:sibson, kwargs...)
    @assert length(x) == length(y) == length(vals) "x, y, and vals must have the same length."
    if !parallel
        for i in eachindex(x, y)
            vals[i] = itp(x[i], y[i], 1; method, kwargs...)
        end
    else
        caches = get_cache(itp)
        nt = length(caches)
        chunked_iterator = chunks(vals, nt)
        Threads.@threads for (xrange, chunk_id) in chunked_iterator
            for i in xrange
                vals[i] = itp(x[i], y[i], chunk_id; method, kwargs...)
            end
        end
    end
    return nothing
end
function (itp::NaturalNeighbourInterpolant)(x::AbstractVector, y::AbstractVector; parallel=true, method=:sibson, kwargs...)
    @assert length(x) == length(y) "x and y must have the same length."
    n = length(x)
    tri = get_triangulation(itp)
    F = number_type(tri)
    vals = zeros(F, n)
    itp(vals, x, y; method, parallel, kwargs...)
    return vals
end

end # module NaturalNeighbourInterp