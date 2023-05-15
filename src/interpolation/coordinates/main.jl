function compute_natural_coordinates(
    method::AbstractInterpolator,
    tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL},
    interpolation_point,
    cache::InterpolantCache{F}=InterpolantCache(tri);
    method=:sibson,
    kwargs...
) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,F}
    if method == Sibson()
        return _compute_sibson_coordinates(tri, interpolation_point, cache; kwargs...)
    elseif method == Triangle()
        return _compute_triangle_coordinates(tri, interpolation_point, cache; kwargs...)
    elseif method == Nearest() # not local coordinates, but still a nice method to have...
        return _compute_nearest_coordinates(tri, interpolation_point, cache; kwargs...)
    elseif method == Laplace()
        return _compute_laplace_coordinates(tri, interpolation_point, cache; kwargs...)
    end
end
