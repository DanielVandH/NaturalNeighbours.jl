abstract type AbstractDifferentiator end
struct Direct <: AbstractDifferentiator end
struct Iterative <: AbstractDifferentiator end
wrap_differentiator(d::AbstractDifferentiator) = d
function wrap_differentiator(d)
    if d == :direct
        return Direct()
    elseif d == :iterative
        return Iterative()
    else
        throw(ArgumentError("Unknown differentiator: $d"))
    end
end

function generate_derivatives(points, z; kwargs...)
    tri = triangulate(points, delete_ghosts=false)
    return generate_derivatives(tri, z; kwargs...)
end
function generate_derivatives(x::AbstractVector, y::AbstractVector, z; kwargs...)
    @assert length(x) == length(y) == length(z) "x, y, and z must have the same length."
    points = [(ξ, η) for (ξ, η) in zip(x, y)]
    return generate_derivatives(points, z; kwargs...)
end
function generate_derivatives(tri::Triangulation, z; kwargs...)
    differentiator = NaturalNeighboursDifferentiator(tri, z)
    return generate_derivatives(differentiator; kwargs...)
end
function generate_derivatives(differentiator::NaturalNeighboursDifferentiator; method=:direct, kwargs...)
    method = wrap_differentiator(method)
    return generate_derivatives(differentiator, method; kwargs...)
end

function generate_derivatives(
    differentiator::NaturalNeighboursDifferentiator,
    method::AbstractDifferentiator;
    parallel=true,
    use_cubic_terms=true,
    use_sibson_weight=false,
    alpha=1 / 2
)
    method = wrap_differentiator(method)
    tri = get_triangulation(differentiator)
    z = get_z(differentiator)
    ∇ = get_gradients(differentiator)
    H = get_hessians(differentiator)
    
end

