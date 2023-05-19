"""
    abstract type AbstractDifferentiator end 

Abstract type for defining the method used for differentiating an interpolant or generating derivatives at data sites. 

- `Direct()`: Generate derivatives directly with one least squares problem.
- `Iterative()`: Generate derivatives iteratively: Gradients are estimated first, and then both gradients and Hessians are estimated with the initial gradients used to refine the results.  
"""
abstract type AbstractDifferentiator end
@doc """
    Direct()

Generate derivatives directly with one least squares problem.
""" struct Direct <: AbstractDifferentiator end
@doc """
    Iterative()

Generate derivatives iteratively: Gradients are estimated first, and then both gradients and Hessians are estimated with the initial gradients used to refine the results.
""" struct Iterative <: AbstractDifferentiator end
dwrap(d::AbstractDifferentiator) = d
function dwrap(d::Symbol)
    if d == :direct
        return Direct()
    elseif d == :iterative
        return Iterative()
    else
        throw(ArgumentError("Unknown differentiator: $d"))
    end
end

"""
    generate_derivatives(
        tri,
        z,
        derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()],
        neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()];
        parallel=true,
        method=Direct(),
        use_cubic_terms=true,
        alpha=0.1,
        initial_gradients=dwrap(method) == Direct() ? nothing : generate_gradients(tri, z, derivative_caches, neighbour_caches; method=dwrap(method), parallel, rng)
    )

Generate derivatives at the data sites defined by the triangulation `tri` with associated function values `tri`.

# Arguments 
- `tri`: A `Triangulation` object.
- `z`: A vector of function values at the data sites.
- `derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()]`: A vector of `DerivativeCache` objects, one for each thread.
- `neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()]`: A vector of `NaturalNeighboursCache` objects, one for each thread.

# Keyword Arguments 
- `parallel=true`: Whether to use multithreading or not.
- `method=Direct()`: The method used for generating the derivatives. See [`AbstractDifferentiator`](@ref).
- `use_cubic_terms=true`: Whether to use cubic terms for estimating the second order derivatives. Only relevant for `method == Direct()`.
- `alpha=0.1`: The weighting parameter used for estimating the second order derivatives. Only relevant for `method == Iterative()`.
- `initial_gradients=dwrap(method) == Direct() ? nothing : generate_gradients(tri, z, derivative_caches, neighbour_caches; method=dwrap(method), parallel, rng)`: The initial gradients used for estimating the second order derivatives. Only relevant for `method == Iterative()`.

# Output 
- `∇`: A vector of gradients at the data sites. Each element is a `Tuple` defining the gradient entries.
- `ℋ`: A vector of Hessians at the data sites. Each element is a `Tuple` defining the Hessian entries in the form `(H[1, 1], H[2, 2], H[1, 2])` (`H[2, 1]` is the same as `H[2, 2]`).
"""
function generate_derivatives(
    tri,
    z,
    derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()],
    neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()];
    parallel=true,
    method=Direct(),
    use_cubic_terms=true,
    alpha=0.1,
    initial_gradients=dwrap(method) == Direct() ? nothing : generate_gradients(tri, z, derivative_caches, neighbour_caches; method=dwrap(method), parallel, rng)
)
    n = length(z)
    F = number_type(tri)
    ∇ = Vector{NTuple{2,F}}(undef, n)
    ℋ = Vector{NTuple{3,F}}(undef, n)
    if !parallel
        generate_second_order_derivatives!(∇, ℋ, method, tri, z, eachindex(z), derivative_caches, neighbour_caches, 1; alpha, use_cubic_terms, initial_gradients)
    else
        nt = length(derivative_caches)
        chunked_iterator = chunks(z, nt)
        Base.Threads.@threads for (zrange, chunk_id) in chunked_iterator
            generate_second_order_derivatives!(∇, ℋ, method, tri, z, zrange, derivative_caches, neighbour_caches, chunk_id; alpha, use_cubic_terms, initial_gradients)
        end
    end
    return ∇, ℋ
end

"""
    generate_gradients(
        tri,
        z,
        derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()],
        neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()];
        parallel=true
    )

Generate gradients at the data sites defined by the triangulation `tri` with associated function values `tri`.

# Arguments
- `tri`: A `Triangulation` object.
- `z`: A vector of function values at the data sites.
- `derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()]`: A vector of `DerivativeCache` objects, one for each thread.
- `neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()]`: A vector of `NaturalNeighboursCache` objects, one for each thread.

# Keyword Arguments
- `parallel=true`: Whether to use multithreading or not.

# Output
- `∇`: A vector of gradients at the data sites. Each element is a `Tuple` defining the gradient entries.
"""
function generate_gradients(
    tri,
    z,
    derivative_caches=[DerivativeCache(tri) for _ in 1:Base.Threads.nthreads()],
    neighbour_caches=[NaturalNeighboursCache(tri) for _ in 1:Base.Threads.nthreads()];
    parallel=true
)
    n = length(z)
    F = number_type(tri)
    ∇ = Vector{NTuple{2,F}}(undef, n)
    if !parallel
        generate_first_order_derivatives!(∇, Direct(), tri, z, eachindex(z), derivative_caches, neighbour_caches, 1)
    else
        nt = length(derivative_caches)
        chunked_iterator = chunks(∇, nt)
        Base.Threads.@threads for (zrange, chunk_id) in chunked_iterator
            generate_first_order_derivatives!(∇, Direct(), tri, z, zrange, derivative_caches, neighbour_caches, chunk_id)
        end
    end
    return ∇
end

# need these redirections to avoid Boxing 

@inline function generate_second_order_derivatives!(∇, ℋ, method, tri, z, zrange, derivative_caches, neighbour_caches, id; alpha, use_cubic_terms, initial_gradients)
    for i in zrange
        zᵢ = z[i]
        d_cache = derivative_caches[id]
        n_cache = neighbour_caches[id]
        S = get_iterated_neighbourhood(d_cache)
        S′ = get_second_iterated_neighbourhood(d_cache)
        if method == Direct()
            λ, E = get_taylor_neighbourhood!(S, S′, tri, i, 2 + use_cubic_terms, n_cache)
        elseif method == Iterative()
            λ, E = get_taylor_neighbourhood!(S, S′, tri, i, 1, n_cache)
        end
        ∇[i], ℋ[i] = generate_second_order_derivatives(method, tri, z, zᵢ, i, λ, E, derivative_caches, id; alpha, use_cubic_terms, initial_gradients)
    end
    return nothing
end

@inline function generate_first_order_derivatives!(∇, method, tri, z, zrange, derivative_caches, neighbour_caches, id)
    for i in zrange
        zᵢ = z[i]
        d_cache = derivative_caches[id]
        n_cache = neighbour_caches[id]
        S = get_iterated_neighbourhood(d_cache)
        S′ = get_second_iterated_neighbourhood(d_cache)
        λ, E = get_taylor_neighbourhood!(S, S′, tri, i, 1, n_cache)
        ∇[i] = generate_first_order_derivatives(method, tri, z, zᵢ, i, λ, E, derivative_caches, id)
    end
    return nothing
end

@inline function generate_first_order_derivatives(
    method::AbstractDifferentiator,
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    derivative_caches,
    id;
    kwargs...)
    return generate_first_order_derivatives(method, tri, z, zᵢ, i, λ, E, derivative_caches[id]; kwargs...)
end
@inline function generate_second_order_derivatives(
    method::AbstractDifferentiator,
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    derivative_caches,
    id;
    kwargs...)
    return generate_second_order_derivatives(method, tri, z, zᵢ, i, λ, E, derivative_caches[id]; kwargs...)
end