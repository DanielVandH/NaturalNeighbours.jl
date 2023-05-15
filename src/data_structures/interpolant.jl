struct NaturalNeighboursolant{T<:Triangulation,F,C}
    triangulation::T
    z::Vector{F}
    cache::C
    function NaturalNeighboursolant(tri::T, z::Vector{F}) where {T,F}
        @assert num_solid_vertices(tri) == length(z) "The number of points in the triangulation must equal the length of the data vector."
        !has_ghost_triangles(tri) && add_ghost_triangles!(tri)
        if has_boundary_nodes(tri)
            throw(ArgumentError("Natural neighbour interpolation is only defined over unconstrained triangulations."))
        end
        nt = Base.Threads.nthreads()
        caches = [InterpolantCache(tri) for _ in 1:nt]
        return new{T,F,typeof(caches)}(tri, z, caches)
    end
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalNeighboursolant)
    z = get_z(nc)
    println(io, "Natural Neighbour Interpolant")
    print(io, "    z: ", z)
end
get_triangulation(ni::NaturalNeighboursolant) = ni.triangulation
get_z(ni::NaturalNeighboursolant) = ni.z
get_cache(ni::NaturalNeighboursolant) = ni.cache
get_cache(ni::NaturalNeighboursolant, id) = ni.cache[id]
