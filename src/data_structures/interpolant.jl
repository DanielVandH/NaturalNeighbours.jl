struct NaturalNeighboursInterpolant{T<:Triangulation,F,G,H,N,D,Z<:Union{<:AbstractVector{F},<:AbstractMatrix{F}}}
    triangulation::T
    z::Z
    gradient::G # (∂ˣf, ∂ʸf)
    hessian::H # (∂ˣˣf, ∂ʸʸf, ∂ˣʸf)
    neighbour_cache::N
    derivative_cache::D
    function NaturalNeighboursInterpolant(
        tri::T,
        z::Z,
        gradient=nothing,
        hessian=nothing;
        derivatives=false,
        kwargs...
    ) where {T,F,Z<:Union{<:AbstractVector{F},<:AbstractMatrix{F}}}
        @assert num_points(tri) == ndata(z) "The number of points in the triangulation must equal the number of function values."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            @warn "Natural neighbour interpolation is only defined over unconstrained triangulations.\nYou may find unexpected results when interpolating near the boundaries or constrained edges, and especially near non-convex boundaries or outside of the triangulation.\nIn your later evaluations of this interpolant, consider using project=false." maxlog = 1
        end
        nt = Base.Threads.nthreads()
        derivative_caches = [DerivativeCache(tri, z) for _ in 1:nt]
        neighbour_caches = [NaturalNeighboursCache(tri) for _ in 1:nt]
        D = typeof(derivative_caches)
        N = typeof(neighbour_caches)
        if derivatives
            ∇, ℋ = generate_derivatives(tri, z, derivative_caches, neighbour_caches; kwargs...)
        else
            ∇ = nothing # TODO: In 2.0, change these to be NTuple{2, F}[] and NTuple{3, F}[]
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
        return new{T,F,G,H,N,D,Z}(tri, z, gradient, hessian, neighbour_caches, derivative_caches)
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
get_z(ni::NaturalNeighboursInterpolant, i) = get_data(get_z(ni), i)
get_neighbour_cache(ni::NaturalNeighboursInterpolant) = ni.neighbour_cache
get_neighbour_cache(ni::NaturalNeighboursInterpolant, id) = get_neighbour_cache(ni)[id]
get_derivative_cache(ni::NaturalNeighboursInterpolant) = ni.derivative_cache
get_derivative_cache(ni::NaturalNeighboursInterpolant, id) = get_derivative_cache(ni)[id]
get_gradient(ni::NaturalNeighboursInterpolant) = ni.gradient
get_gradient(ni::NaturalNeighboursInterpolant, i) = get_data(get_gradient(ni), i)
get_hessian(ni::NaturalNeighboursInterpolant) = ni.hessian
get_hessian(ni::NaturalNeighboursInterpolant, i) = get_data(get_hessian(ni), i)
has_gradient(ni::NaturalNeighboursInterpolant) = !isnothing(get_gradient(ni))
has_hessian(ni::NaturalNeighboursInterpolant) = !isnothing(get_hessian(ni))