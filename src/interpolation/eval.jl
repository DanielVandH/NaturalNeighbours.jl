function _eval_natural_coordinates(coordinates, indices, z)
    val = zero(eltype(z))
    for (λ, k) in zip(coordinates, indices)
        zₖ = z[k]
        val += λ * zₖ
    end
    return val
end

function _eval_natural_coordinates(nc::NaturalCoordinates{F}, z) where {F}
    coordinates = get_coordinates(nc)
    indices = get_indices(nc)
    return _eval_natural_coordinates(coordinates, indices, z)
end

function _eval_interp(method, itp::NaturalNeighboursInterpolant, p, cache; kwargs...)
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(method, tri, p, cache; kwargs...)
    z = get_z(itp)
    return _eval_natural_coordinates(nc, z)
end

function _eval_natural_coordinates(::Sibson{1}, nc::NaturalCoordinates{F}, z, gradients, tri) where {F}
    sib0 = _eval_natural_coordinates(nc, z)
    coordinates = get_coordinates(nc)
    if length(coordinates) ≤ 2 # 2 means extrapolation, 1 means we're evaluating at a data site 
        return sib0
    end
    ζ, α, β = _compute_sibson_1_coordinates(nc, tri, z, gradients)
    num = α * sib0 + β * ζ
    den = α + β
    return num / den
end

function _eval_interp(method::Sibson{1}, itp::NaturalNeighboursInterpolant, p, cache; kwargs...) # # has to be a different form since Sib0 blends two functions 
    gradients = get_gradient(itp)
    if isnothing(gradients)
        throw(ArgumentError("Gradients must be provided for Sibson-1 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true)."))
    end
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(Sibson(), tri, p, cache; kwargs...)
    z = get_z(itp)
    return _eval_natural_coordinates(method, nc, z, gradients, tri)
end

