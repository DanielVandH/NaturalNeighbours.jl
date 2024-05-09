function _compute_nearest_coordinates(
    tri::Triangulation{P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE},
    interpolation_point,
    cache::NaturalNeighboursCache{F}=NaturalNeighboursCache(tri);
    kwargs...
) where {P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE,F}
    coordinates = get_coordinates(cache)
    envelope = get_envelope(cache)
    last_triangle = get_last_triangle(cache)
    i = jump_to_voronoi_polygon(tri, interpolation_point; try_points=last_triangle[])
    resize!(coordinates, 1)
    resize!(envelope, 1)
    coordinates[1] = one(number_type(tri))
    envelope[1] = i
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end

function compute_natural_coordinates(::Nearest, tri, interpolation_point, cache=NaturalNeighboursCache(tri); kwargs...)
    return _compute_nearest_coordinates(tri, interpolation_point, cache; kwargs...)
end