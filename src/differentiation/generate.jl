function generate_derivatives(points, z, order; parallel = true, kwargs...)
    tri = triangulate(points, delete_ghosts=false; kwargs...)
    return generate_derivatives(tri, z, order; parallel, kwargs...)
end
function generate_derivatives(x::AbstractVector, y::AbstractVector, z, order; method=:direct, parallel=true, kwargs...)
    @assert length(x) == length(y) == length(z) "x, y, and z must have the same length."
    points = [(ξ, η) for (ξ, η) in zip(x, y)]
    return generate_derivatives(points, z, order; parallel, kwargs...)
end
function generate_derivatives(tri::Triangulation, z, order; parallel=true, kwargs...)
    differentiator = NaturalNeighboursDifferentiator(tri, z)
    return generate_derivatives(differentiator, z; parallel, kwargs...)
end

function generate_derivatives(differentiator::NaturalNeighboursDifferentiator, order; method=:direct, parallel=true, kwargs...)
    @assert order ∈ (1, 2) "Only gradients and Hessians can be estimated."
    if order == 1
        return generate_gradients(differentiator; parallel, kwargs...)
    else
        return generate_gradients_and_hessians(differentiator; parallel, kwargs...)
    end
end

function generate_gradients(differentiator; method=:direct, parallel=true, kwargs...)
    tri = get_triangulation(differentiator)
    z = get_z(differentiator)
    caches = get_cache(differentiator)
    F = number_type(tri)
    gradients = Vector{NTuple{3, F}}(undef, length(z))
    if !parallel
        cache = get_cache(differentiator, 1)
        for i in eachindex(z)
            p = get_point(tri, i)
            x, y = getxy(p)
            gradients[i] = _eval_gradient(differentiator, p, cache; kwargs...)

end