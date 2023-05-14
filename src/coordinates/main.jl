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
    elseif method == :nearest # not local coordinates, but still a nice method to have...
        return _compute_nearest_coordinates(tri, interpolation_point, cache; kwargs...)
    elseif method == :laplace
        return _compute_laplace_coordinates(tri, interpolation_point, cache; kwargs...)
    else
        throw(ArgumentError("method must be one of :sibson, :triangle, :laplace, or :nearest."))
    end
end
