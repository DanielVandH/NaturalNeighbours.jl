struct NaturalNeighboursDifferentiator{T,F,C}
    triangulation::T 
    z::Vector{F}
    cache::C
    function NaturalNeighboursDifferentiator(tri::T, z::Vector{F}) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        caches = [DerivativeCache(tri) for _ in 1:nt]
        return new{T,F,typeof(caches)}(tri, z, caches)
    end
end
get_triangulation(nd::NaturalNeighboursDifferentiator) = nd.triangulation
get_z(nd::NaturalNeighboursDifferentiator) = nd.z
get_z(nd::NaturalNeighboursDifferentiator, i) = nd.z[i]
get_cache(nd::NaturalNeighboursDifferentiator) = nd.cache
get_cache(nd::NaturalNeighboursDifferentiator, id) = nd.cache[id]