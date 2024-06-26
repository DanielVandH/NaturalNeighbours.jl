function _compute_laplace_coordinates(
    tri::Triangulation{P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE},
    interpolation_point,
    cache::NaturalNeighboursCache{F}=NaturalNeighboursCache(tri);
    project=true,
    kwargs...
) where {P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE,F}
    coordinates = get_coordinates(cache)
    envelope = get_envelope(cache)
    insertion_event_history = get_insertion_event_history(cache)
    temp_adjacent = get_temp_adjacent(cache)
    last_triangle = get_last_triangle(cache)
    envelope, temp_adjacent, insertion_event_history, V = compute_bowyer_envelope!(envelope, tri, insertion_event_history, temp_adjacent, interpolation_point; try_points=last_triangle[], kwargs...) #kwargs are add_point! kwargs
    i, j, return_flag = check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    return_flag && return two_point_interpolate!(coordinates, envelope, tri, i, j, interpolation_point, project)
    resize!(coordinates, length(envelope) - 1)
    w = zero(number_type(tri))
    for i in firstindex(envelope):(lastindex(envelope)-1)
        ratio, u, prev_u, next_u = laplace_ratio(tri, envelope, i, interpolation_point) # could reuse a circumcenter here, but it's not the dominating part of the computation anyway.
        isnan(ratio) && return handle_duplicate_points!(tri, interpolation_point, coordinates, envelope, u, prev_u, next_u)
        coordinates[i] = max(ratio, zero(ratio)) # coordinate types like Float32 can sometimes get e.g. -1f-8
        w += coordinates[i]
    end
    pop!(envelope)
    coordinates ./= w
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end

function laplace_ratio(tri, envelope, i, interpolation_point)
    u = envelope[i]
    prev_u = envelope[previndex_circular(envelope, i)]
    next_u = envelope[nextindex_circular(envelope, i)]
    p, q, r = get_point(tri, u, prev_u, next_u)
    g1 = triangle_circumcenter(q, p, interpolation_point)
    g2 = triangle_circumcenter(p, r, interpolation_point)
    g1x, g1y = getxy(g1)
    g2x, g2y = getxy(g2)
    ℓ² = (g1x - g2x)^2 + (g1y - g2y)^2
    px, py = getxy(p)
    x, y = getxy(interpolation_point)
    d² = (px - x)^2 + (py - y)^2
    w = sqrt(ℓ² / d²)
    return w, u, prev_u, next_u
end

function compute_natural_coordinates(::Laplace, tri, interpolation_point, cache=NaturalNeighboursCache(tri); kwargs...) 
    return _compute_laplace_coordinates(tri, interpolation_point, cache; kwargs...)
end