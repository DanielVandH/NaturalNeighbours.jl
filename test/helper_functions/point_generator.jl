function random_points_in_convex_hull(tri::Triangulation, n; rng=Random.default_rng()) # bit slow. oh well
    boundary_nodes = get_convex_hull_vertices(tri)
    points = get_points(tri)
    bbox = DT.polygon_bounds(points, boundary_nodes)
    F = DT.number_type(tri)
    pts = NTuple{2,F}[]
    while length(pts) < n
        p = (rand(rng, F) * (bbox[2] - bbox[1]) + bbox[1], rand(rng, F) * (bbox[4] - bbox[3]) + bbox[3])
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        if δ > 0
            push!(pts, p)
        end
    end
    return pts
end