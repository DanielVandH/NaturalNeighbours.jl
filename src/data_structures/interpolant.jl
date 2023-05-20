struct NaturalNeighboursInterpolant{T<:Triangulation,F,G,H,N,D}
    triangulation::T
    z::Vector{F}
    gradient::G # (∂ˣf, ∂ʸf)
    hessian::H # (∂ˣˣf, ∂ʸʸf, ∂ˣʸf)
    neighbour_cache::N
    derivative_cache::D
    function NaturalNeighboursInterpolant(
        tri::T,
        z::AbstractVector{F},
        gradient=nothing,
        hessian=nothing;
        derivatives=false,
        kwargs...
    ) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        derivative_caches = [DerivativeCache(tri) for _ in 1:nt]
        neighbour_caches = [NaturalNeighboursCache(tri) for _ in 1:nt]
        D = typeof(derivative_caches)
        N = typeof(neighbour_caches)
        if derivatives
            ∇, ℋ = generate_derivatives(tri, z, derivative_caches, neighbour_caches; kwargs...)
        else
            ∇ = nothing
            ℋ = nothing
        end
        if isnothing(gradient)
            gradient = ∇
        end
        if isnothing(hessian)
            hessian = ℋ
        end
        G = typeof(gradient)
        H = typeof(hessian)
        return new{T,F,G,H,N,D}(tri, z, gradient, hessian, neighbour_caches, derivative_caches)
    end
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalNeighboursInterpolant)
    z = get_z(nc)
    println(io, "Natural Neighbour Interpolant")
    println(io, "    z: ", z)
    println(io, "    ∇: ", get_gradient(nc))
    print(io, "    H: ", get_hessian(nc))
end
get_triangulation(ni::NaturalNeighboursInterpolant) = ni.triangulation
get_z(ni::NaturalNeighboursInterpolant) = ni.z
get_z(ni::NaturalNeighboursInterpolant, i) = ni.z[i]
get_neighbour_cache(ni::NaturalNeighboursInterpolant) = ni.neighbour_cache
get_neighbour_cache(ni::NaturalNeighboursInterpolant, id) = ni.neighbour_cache[id]
get_derivative_cache(ni::NaturalNeighboursInterpolant) = ni.derivative_cache
get_derivative_cache(ni::NaturalNeighboursInterpolant, id) = ni.derivative_cache[id]
get_gradient(ni::NaturalNeighboursInterpolant) = ni.gradient
get_gradient(ni::NaturalNeighboursInterpolant, i) = ni.gradient[i]
get_hessian(ni::NaturalNeighboursInterpolant) = ni.hessian
get_hessian(ni::NaturalNeighboursInterpolant, i) = ni.hessian[i]
has_gradient(ni::NaturalNeighboursInterpolant) = !isnothing(get_gradient(ni))
has_hessian(ni::NaturalNeighboursInterpolant) = !isnothing(get_hessian(ni))