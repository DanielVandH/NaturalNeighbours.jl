interpolate(tri::Triangulation, z; kwargs...) = NaturalNeighbourInterpolant(tri, z)
function interpolate(points, z; kwargs...)
    tri = triangulate(points, delete_ghosts=false; kwargs...)
    return interpolate(tri, z)
end
function interpolate(x::AbstractVector, y::AbstractVector, z; kwargs...)
    @assert length(x) == length(y) == length(z) "x, y, and z must have the same length."
    points = [(ξ, η) for (ξ, η) in zip(x, y)]
    return interpolate(points, z; kwargs...)
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
