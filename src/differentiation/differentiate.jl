"""
    differentiate(itp::NaturalNeighboursInterpolant, order)

Differentiate a given interpolant `itp` up to degree `order` (1 or 2). The returned object is a 
`NaturalNeighboursDifferentiator` struct, which is callable. 

For calling the resulting struct, we define the following methods:

    (∂::NaturalNeighboursDifferentiator)(x, y, zᵢ, nc, id::Integer=1; parallel=false, method=default_diff_method(∂), kwargs...)
    (∂::NaturalNeighboursDifferentiator)(x, y, id::Integer=1; parallel=false, method=default_diff_method(∂), interpolant_method=Sibson(), rng=Random.default_rng(), kwargs...)
    (∂::NaturalNeighboursDifferentiator)(vals::AbstractVector, x::AbstractVector, y::AbstractVector; parallel=true, method=default_diff_method(∂), interpolant_method=Sibson(), kwargs...)
    (∂::NaturalNeighboursDifferentiator{I, O})(x::AbstractVector, y::AbstractVector; parallel=true, method=default_diff_method(∂), interpolant_method=Sibson(), kwargs...) where {I, O}

1. This method is useful if you already have an estimate for the function value, `zᵢ`, at the data site, `(x, y)`, provided you also provide the `NaturalCoordinates` `nc`. `id` is the thread id.
2. This method is for scalars, with `id` referring to a thread id.
3. This method is an in-place method for vectors, storing `∂(x[i], y[i], 1)` into `vals[i]`.
4. This method is similar to (3), but `vals` is constructed and returned.

The available keyword arguments are:
- `parallel=true`: Whether to use multithreading. Ignored for the first two methods. 
- `method=default_diff_method(∂)`: Default method for evaluating the interpolant. `default_diff_method(∂)` returns `Direct()` if the underlying interpolant has no gradients, and `Iterative()` otherwise. The method must be a [`AbstractDifferentiator`](@ref).
- `interpolant_method=Sibson()`: The method used for evaluating the interpolant to estimate `zᵢ` for the latter three methods. See [`AbstractInterpolator`](@ref) for the avaiable methods.
- `rng=Random.default_rng()`: The random number generator used for estimating `zᵢ` for the latter three methods, or for constructing the natural coordinates.
- `use_cubic_terms=true`: If estimating second order derivatives, whether to use cubic terms. Only relevant for `method == Direct()`.
- `alpha=0.1`: If estimating second order derivatives, the weighting parameter used for estimating the second order derivatives. Only relevant for `method == Iterative()`.
- `use_sibson_weight=true`: Whether to weight the residuals in the associated least squares problems by the associated Sibson coordinates. Only relevant for `method == Iterative()` if `order == 2`.

The outputs are:
- `order == 1`: The scalar methods return a `Tuple` of the form `(∂x, ∂y)`, while the vector methods return a vector of `Tuple`s of the form `(∂x, ∂y)`.
- `order == 2`: The scalar methods return a `(∇, ℋ)`, where `∇` is a `Tuple` of the form `(∂x, ∂y)` and `ℋ` is a `Tuple` of the form `(∂xx, ∂yy, ∂xy)`. The vector methods return a vector of `(∇, ℋ)`s.
"""
differentiate(itp::NaturalNeighboursInterpolant, order) = NaturalNeighboursDifferentiator(itp, order)

function _eval_differentiator(method::AbstractDifferentiator, ∂::NaturalNeighboursDifferentiator{I,O}, p, zᵢ, nc, id;
    use_cubic_terms=true,
    alpha=0.1,
    use_sibson_weight=true) where {I,O}
    itp = get_interpolant(∂)
    tri = get_triangulation(itp)
    z = get_z(itp)
    d_cache = get_derivative_cache(itp, id)
    initial_gradients = get_gradient(itp)
    S = get_iterated_neighbourhood(d_cache)
    S′ = get_second_iterated_neighbourhood(d_cache)
    if O == 1
        λ, E = get_taylor_neighbourhood!(S, S′, tri, 1, nc)
        return generate_first_order_derivatives(method, tri, z, zᵢ, p, λ, E, d_cache; use_cubic_terms, alpha, use_sibson_weight, initial_gradients)
    else # O == 2
        λ, E = get_taylor_neighbourhood!(S, S′, tri, 2 + use_cubic_terms, nc)
        return generate_second_order_derivatives(method, tri, z, zᵢ, p, λ, E, d_cache; use_cubic_terms, alpha, use_sibson_weight, initial_gradients)
    end
end

default_diff_method(∂) = isnothing(get_gradient(get_interpolant(∂))) ? Direct() : Iterative()

function (∂::NaturalNeighboursDifferentiator)(x, y, zᵢ, nc, id::Integer=1; parallel=false, method=default_diff_method(∂), kwargs...)
    method = dwrap(method)
    p = (x, y)
    return _eval_differentiator(method, ∂, p, zᵢ, nc, id; kwargs...)
end
function (∂::NaturalNeighboursDifferentiator)(x, y, id::Integer=1; parallel=false, method=default_diff_method(∂), interpolant_method=Sibson(), rng=Random.default_rng(), kwargs...)
    method = dwrap(method)
    interpolant_method = iwrap(interpolant_method)
    p = (x, y)
    itp = get_interpolant(∂)
    tri = get_triangulation(itp)
    n_cache = get_neighbour_cache(itp, id)
    z = get_z(itp)
    if interpolant_method == Sibson(1)
        nc = compute_natural_coordinates(Sibson(), tri, p, n_cache; rng)
        gradients = get_gradient(itp)
        zᵢ = _eval_natural_coordinates(nc, z, gradients, tri)
    else
        nc = compute_natural_coordinates(interpolant_method, tri, p, n_cache; rng)
        zᵢ = _eval_natural_coordinates(nc, z)
    end
    if interpolant_method == Triangle() || interpolant_method == Nearest() # coordinates need to be the natural neighbours
        nc = compute_natural_coordinates(Sibson(), tri, p, n_cache; rng)
    end
    return ∂(x, y, zᵢ, nc, id; parallel, method, kwargs...)
end
function (∂::NaturalNeighboursDifferentiator)(vals::AbstractVector, x::AbstractVector, y::AbstractVector; parallel=true, method=default_diff_method(∂), interpolant_method=Sibson(), kwargs...)
    @assert length(x) == length(y) == length(vals) "x, y, and vals must have the same length."
    method = dwrap(method)
    interpolant_method = iwrap(interpolant_method)
    if !parallel
        for i in eachindex(x, y)
            vals[i] = ∂(x[i], y[i], 1; method, interpolant_method, kwargs...)
        end
    else
        nt = Base.Threads.nthreads()
        chunked_iterator = chunks(vals, nt)
        Threads.@threads for (xrange, chunk_id) in chunked_iterator
            for i in xrange
                vals[i] = ∂(x[i], y[i], chunk_id; method, interpolant_method, kwargs...)
            end
        end
    end
    return nothing
end
function (∂::NaturalNeighboursDifferentiator{I,O})(x::AbstractVector, y::AbstractVector; parallel=true, method=default_diff_method(∂), interpolant_method=Sibson(), kwargs...) where {I,O}
    @assert length(x) == length(y) "x and y must have the same length."
    n = length(x)
    itp = get_interpolant(∂)
    tri = get_triangulation(itp)
    F = number_type(tri)
    if O == 1
        vals = Vector{NTuple{2,F}}(undef, n)
    else # O == 2
        vals = Vector{Tuple{NTuple{2,F},NTuple{3,F}}}(undef, n)
    end
    method = dwrap(method)
    interpolant_method = iwrap(interpolant_method)
    ∂(vals, x, y; method, interpolant_method, parallel, kwargs...)
    return vals
end