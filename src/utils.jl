"""
    identify_exterior_points(x, y, points, boundary_nodes; tol = 0.0)

Given a polygon described by `(points, boundary_nodes)`, matching the 
specification of polygons in DelaunayTriangulation.jl (see [here](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/boundary_handling/)),
returns a vector of indices of the points defined by `(x, y)` that are outside of the polygon.
"""
function identify_exterior_points(x, y, points, boundary_nodes; tol = 0.0)
    @assert length(x) == length(y) "x and y must have the same length."
    exterior_points = Int64[]
    sizehint!(exterior_points, isqrt(length(x)))
    for i in eachindex(x, y)
        xᵢ = x[i]
        yᵢ = y[i]
        q = (xᵢ, yᵢ)
        δ = distance_to_polygon(q, points, boundary_nodes)
        if δ < tol
            push!(exterior_points, i)
        end
    end
    return exterior_points
end