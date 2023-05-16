struct NaturalNeighboursDifferentiator{T,F,C,N}
    triangulation::T
    z::Vector{F}
    gradients::Vector{NTuple{2, F}}
    hessians::Vector{NTuple{3, F}}
    derivative_cache::C
    neighbour_cache::N
    function NaturalNeighboursDifferentiator(tri::T, z::Vector{F}) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        derivative_caches = [DerivativeCache(tri) for _ in 1:nt]
        neighbour_caches = [NaturalNeighboursCache(tri) for _ in 1:nt]
        C = typeof(derivative_caches)
        N = typeof(neighbour_caches)
        n = length(z)
        gradients = Vector{NTuple{2, F}}(undef, n)
        hessians = Vector{NTuple{3, F}}(undef, n)
        return new{T,F,C,N}(tri, z, gradients, hessians, derivative_caches, neighbour_caches)
    end
end
get_triangulation(nd::NaturalNeighboursDifferentiator) = nd.triangulation
get_z(nd::NaturalNeighboursDifferentiator) = nd.z
get_z(nd::NaturalNeighboursDifferentiator, i) = nd.z[i]
get_gradients(nd::NaturalNeighboursDifferentiator) = nd.gradients
get_gradients(nd::NaturalNeighboursDifferentiator, i) = nd.gradients[i]
get_hessians(nd::NaturalNeighboursDifferentiator) = nd.hessians
get_hessians(nd::NaturalNeighboursDifferentiator, i) = nd.hessians[i]
get_derivative_cache(nd::NaturalNeighboursDifferentiator) = nd.derivative_cache
get_derivative_cache(nd::NaturalNeighboursDifferentiator, id) = get_derivative_cache(nd)[id]
get_neighbour_cache(nd::NaturalNeighboursDifferentiator) = nd.neighbour_cache
get_neighbour_cache(nd::NaturalNeighboursDifferentiator, id) = get_neighbour_cache(nd)[id]
