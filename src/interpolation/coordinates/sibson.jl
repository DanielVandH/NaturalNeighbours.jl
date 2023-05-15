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
    for i in firstindex(envelope):(lastindex(envelope)-1)
        pre = pre_insertion_area!(poly_points, envelope, i, tri)
        post = post_insertion_area(envelope, i, tri, interpolation_point)
        isnan(post) && return handle_duplicate_points!(tri, interpolation_point, coordinates, envelope)
        coordinates[i] = pre - post
        w += coordinates[i]
    end
    pop!(envelope)
    coordinates ./= w
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
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

function compute_natural_coordinates(::Sibson, tri, interpolation_point, cache=InterpolantCache(tri); kwargs...) 
    return _compute_sibson_coordinates(tri, interpolation_point, cache; kwargs...)
end