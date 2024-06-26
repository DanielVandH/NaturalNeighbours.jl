@inline function _eval_natural_coordinates(coordinates, indices, z)
    val = zero(eltype(z))
    for (λ, k) in zip(coordinates, indices)
        zₖ = z[k]
        val += λ * zₖ
    end
    return val
end

@inline function _eval_natural_coordinates(nc::NaturalCoordinates{F}, z) where {F}
    coordinates = get_coordinates(nc)
    indices = get_indices(nc)
    return _eval_natural_coordinates(coordinates, indices, z)
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

@inline function _eval_natural_coordinates(::Farin, nc::NaturalCoordinates{F}, z, gradients, tri) where {F}
    λ = get_coordinates(nc)
    if length(λ) ≤ 2 # 2 means extrapolation, 1 means we're evaluating at a data site 
        return _eval_natural_coordinates(nc, z)
    end
    return _compute_farin_coordinates(nc, tri, z, gradients)
end

@inline function _eval_natural_coordinates(::Hiyoshi{2}, nc::NaturalCoordinates{F}, z, gradients, hessians, tri) where {F}
    return _compute_hiyoshi_coordinates(nc, tri, z, gradients, hessians)
end

@inline function _eval_interp(method, itp::NaturalNeighboursInterpolant, p, cache; kwargs...)
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(method, tri, p, cache; kwargs...)
    z = get_z(itp)
    return _eval_natural_coordinates(nc, z)
end

@inline function _eval_interp(method::Union{<:Farin,Sibson{1}}, itp::NaturalNeighboursInterpolant, p, cache; kwargs...)
    gradients = get_gradient(itp)
    if isnothing(gradients)
        throw(ArgumentError("Gradients must be provided for Sibson-1, Farin, or Hiyoshi-2 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true)."))
    end
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(Sibson(), tri, p, cache; kwargs...)
    z = get_z(itp)
    λ = get_coordinates(nc)
    if length(λ) ≤ 2 # 2 means extrapolation, 1 means we're evaluating at a data site 
        return _eval_natural_coordinates(nc, z)
    end
    return _eval_natural_coordinates(method, nc, z, gradients, tri)
end

@inline function _eval_interp(method::Hiyoshi{2}, itp::NaturalNeighboursInterpolant, p, cache; kwargs...)
    gradients = get_gradient(itp)
    hessians = get_hessian(itp)
    if isnothing(gradients) || isnothing(hessians)
        throw(ArgumentError("Gradients and Hessians must be provided for Hiyoshi-2 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true)."))
    end
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(Sibson(), tri, p, cache; kwargs...)
    z = get_z(itp)
    λ = get_coordinates(nc)
    if length(λ) ≤ 2 # 2 means extrapolation, 1 means we're evaluating at a data site 
        return _eval_natural_coordinates(nc, z)
    end
    return _eval_natural_coordinates(method, nc, z, gradients, hessians, tri)
end

function _get_nc_and_z(method::AbstractInterpolator{D}, p, z, gradients, hessians, tri, cache=NaturalNeighboursCache(tri); rng=Random.default_rng(), project=true) where {D}
    if method == Triangle() || method == Nearest() # coordinates need to be the natural neighbours
        nc = compute_natural_coordinates(Sibson(), tri, p, cache; rng, project)
    else
        nc = compute_natural_coordinates(method, tri, p, cache; rng, project)
    end
    if D == 0
        zᵢ = _eval_natural_coordinates(nc, z)
    elseif D == 1
        zᵢ = _eval_natural_coordinates(method, nc, z, gradients, tri)
    else # D == 2 
        zᵢ = _eval_natural_coordinates(method, nc, z, gradients, hessians, tri)
    end
    return nc, zᵢ
end