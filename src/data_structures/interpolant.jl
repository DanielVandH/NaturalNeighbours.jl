struct NaturalNeighboursInterpolant{T<:Triangulation,F,C,G,H}
    triangulation::T
    z::Vector{F}
    gradient::G # (∂ˣf, ∂ʸf)
    hessian::H # (∂ˣˣf, ∂ʸʸf, ∂ˣʸf)
    cache::C
    function NaturalNeighboursInterpolant(tri::T, z::Vector{F}, gradient=nothing, hessian=nothing) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        caches = [InterpolantCache(tri) for _ in 1:nt]
        return new{T,F,typeof(caches),typeof(gradient),typeof(hessian)}(tri, z, gradient, hessian, caches)
    end
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalNeighboursInterpolant)
    z = get_z(nc)
    println(io, "Natural Neighbour Interpolant")
    print(io, "    z: ", z)
end
get_triangulation(ni::NaturalNeighboursInterpolant) = ni.triangulation
get_z(ni::NaturalNeighboursInterpolant) = ni.z
get_z(ni::NaturalNeighboursInterpolant, i) = ni.z[i]
get_cache(ni::NaturalNeighboursInterpolant) = ni.cache
get_cache(ni::NaturalNeighboursInterpolant, id) = ni.cache[id]
get_gradient(ni::NaturalNeighboursInterpolant) = ni.gradient
get_gradient(ni::NaturalNeighboursInterpolant, i) = ni.gradient[i]
get_hessian(ni::NaturalNeighboursInterpolant) = ni.hessian
get_hessian(ni::NaturalNeighboursInterpolant, i) = ni.hessian[i]
has_gradient(ni::NaturalNeighboursInterpolant) = !isnothing(get_gradient(ni))
has_hessian(ni::NaturalNeighboursInterpolant) = !isnothing(get_hessian(ni))