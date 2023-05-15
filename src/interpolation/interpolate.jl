abstract type AbstractInterpolator end
struct Sibson <: AbstractInterpolator end
struct Triangle <: AbstractInterpolator end
struct Nearest <: AbstractInterpolator end
struct Laplace <: AbstractInterpolator end
wrap_interpolator(s::AbstractInterpolator) = s
function wrap_interpolator(s)
    if s == :sibson
        return Sibson()
    elseif s == :triangle
        return Triangle()
    elseif s == :nearest
        return Nearest()
    elseif s == :laplace
        return Laplace()
    else
        throw(ArgumentError("Unknown interpolator: $s"))
    end
end

interpolate(tri::Triangulation, z; kwargs...) = NaturalNeighboursInterpolant(tri, z)
function interpolate(points, z; kwargs...)
    tri = triangulate(points, delete_ghosts=false; kwargs...)
    return interpolate(tri, z)
end
function interpolate(x::AbstractVector, y::AbstractVector, z; kwargs...)
    @assert length(x) == length(y) == length(z) "x, y, and z must have the same length."
    points = [(ξ, η) for (ξ, η) in zip(x, y)]
    return interpolate(points, z; kwargs...)
end

function _eval_interp(method, itp::NaturalNeighboursInterpolant, p, cache; kwargs...)
    tri = get_triangulation(itp)
    nc = compute_natural_coordinates(method, tri, p, cache; kwargs...)
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

function (itp::NaturalNeighboursInterpolant)(x, y, id::Integer=1; parallel=false, method=Sibson(), kwargs...)
    p = (x, y)
    cache = get_cache(itp, id)
    return _eval_interp(wrap_interpolator(method), itp, p, cache; kwargs...)
end

function (itp::NaturalNeighboursInterpolant)(vals::AbstractVector, x::AbstractVector, y::AbstractVector; parallel=true, method=Sibson(), kwargs...)
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
function (itp::NaturalNeighboursInterpolant)(x::AbstractVector, y::AbstractVector; parallel=true, method=Sibson(), kwargs...)
    @assert length(x) == length(y) "x and y must have the same length."
    n = length(x)
    tri = get_triangulation(itp)
    F = number_type(tri)
    vals = zeros(F, n)
    itp(vals, x, y; method, parallel, kwargs...)
    return vals
end
