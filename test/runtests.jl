using NaturalNeighbourInterp
using Test

@testset "NaturalNeighbourInterp.jl" begin
    # Write your tests here.
end

using DelaunayTriangulation, Random, CairoMakie, StableRNGs
const DT = DelaunayTriangulation
const NNI = NaturalNeighbourInterp

function random_points_in_convex_hull(tri::Triangulation, n) # bit slow. oh well
    boundary_nodes = get_convex_hull_indices(tri)
    points = get_points(tri)
    bbox = DT.polygon_bounds(points, boundary_nodes)
    F = DT.number_type(tri)
    pts = NTuple{2,F}[]
    while length(pts) < n
        p = (rand(F) * (bbox[2] - bbox[1]) + bbox[1], rand(F) * (bbox[4] - bbox[3]) + bbox[3])
        δ = DT.distance_to_polygon(p, points, boundary_nodes)
        if δ > 0
            push!(pts, p)
        end
    end
    return pts
end

@testset "Computing the Bowyer-Watson envelope" begin
    pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
    tri = triangulate(pts, randomise=false, delete_ghosts=false)
    envelope = Int64[]
    history = DT.initialise_event_history(tri)
    A = DT.Adjacent{Int64,NTuple{2,Int64}}()
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (6.0, 4.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 4.0)))
    @test DT.circular_equality(envelope, [1, 5, 7, 6, 1])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (6.0, 1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 1.0)))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (7.0, -1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (7.0, -1.0)))
    @test DT.circular_equality(envelope, [2, DT.BoundaryIndex, 3, 8, 7, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (0.0, 3.0))
    @test DT.is_on(DT.point_position_relative_to_triangle(tri, V, (0.0, 3.0)))
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (6.0, 4.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 4.0)))
    for T in DT.each_added_triangle(history)
        i, j, k = indices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A) == 3length(DT.each_added_triangle(history))
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (6.0, 1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 1.0)))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    for T in DT.each_added_triangle(history)
        i, j, k = indices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A) == 3length(DT.each_added_triangle(history))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (7.0, -1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (7.0, -1.0)))
    for T in DT.each_added_triangle(history)
        i, j, k = indices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A) == 3length(DT.each_added_triangle(history))
    @test DT.circular_equality(envelope, [2, DT.BoundaryIndex, 3, 8, 7, 2])
end

@testset "polygon_area" begin
    θ = LinRange(0, 2π - 0.1, 25)
    z = exp.(im * θ)
    points = [(real(z[i]), imag(z[i])) for i in eachindex(z)]
    push!(points, points[begin])
    A = NNI.polygon_area(points)
    _A, _ = DT.polygon_features(points, [1:length(points)-1; 1])
    @test A ≈ _A

    pts = rand(2, 50)
    tri = triangulate(pts)
    boundary_nodes = get_convex_hull_indices(tri)
    A = NNI.polygon_area(pts[:, boundary_nodes])
    _A = DT.get_total_area(tri)
    @test A ≈ _A
end

@testset "Natural coordinates for interior points" begin
    pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
    tri = triangulate(pts, randomise=false, delete_ghosts=false)
    n = 2500
    pts = random_points_in_convex_hull(tri, n)
    for p in pts
        natural_coordinates = NNI.compute_natural_coordinates(tri, p)
        δ = NNI.get_barycentric_deviation(natural_coordinates)
        @test δ ≈ 0 atol = 1e-9
    end

    for _ in 1:50
        pts = [(randn() + rand(), rand() - 2randn()) for _ in 1:250]
        tri = triangulate(pts, delete_ghosts=false)
        n = 5000
        random_points = random_points_in_convex_hull(tri, n)
        for p in random_points
            natural_coordinates = NNI.compute_natural_coordinates(tri, p)
            δ = NNI.get_barycentric_deviation(natural_coordinates)
            @test δ ≈ 0 atol = 1e-5
        end
    end
end

@testset "Two-point interpolations" begin
    for _ in 1:10
        tri = triangulate(rand(2, 500))
        coordinates = zeros(5)
        envelope = zeros(Int, 5)
        for _ in 1:100
            e = (rand ∘ get_edges)(tri)
            i, j = DT.edge_indices(e)
            t = rand()
            p = (1 - t) .* get_point(tri, i) .+ t .* get_point(tri, j)
            nc = NNI.two_point_interpolate!(coordinates, envelope, tri, i, j, p)
            @test sum(NNI.get_coordinates(nc)) ≈ 1
            @test NNI.get_indices(nc) == [i, j]
            @test NNI.get_interpolation_point(nc) == p
            @test NNI.get_triangulation(nc) == tri
            @test NNI.get_coordinates(nc) ≈ [1 - t, t]
            λ, k = NNI.get_coordinates(nc), NNI.get_indices(nc)
            @test collect(p) ≈ collect(λ[1] .* get_point(tri, k[1]) .+ λ[2] .* get_point(tri, k[2]))
            @test NNI.get_barycentric_deviation(nc) ≈ 0 atol = 1e-8
        end
    end
end

@testset "Interpolation" begin
    rng = StableRNG(123)
    pts = [(rand(rng), rand(rng)) for _ in 1:50]
    tri = triangulate(pts, rng=rng, delete_ghosts=false)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    z = [f(x, y) for (x, y) in pts]

    itp = interpolate(tri, z)
    @test DT.get_triangulation(itp) == tri
    @test NNI.get_z(itp) == z
    @test length(NNI.get_cache(itp)) == Base.Threads.nthreads()
    @test NNI.get_cache(itp, 1) == itp.cache[1]
    @test itp isa NNI.NaturalNeighbourInterpolant
    DT.lock_convex_hull!(tri)
    @test_throws ArgumentError interpolate(tri, z)
    DT.unlock_convex_hull!(tri)
    @test_throws AssertionError interpolate(tri, z[1:end-1])

    x = getx.(pts)
    y = gety.(pts)
    for _ in 1:500
        vals = itp(x, y; parallel=false)
        vals2 = similar(vals)
        itp(vals2, x, y; parallel=false)
        vals3 = itp(x, y; parallel=true)
        vals4 = similar(vals3)
        itp(vals4, x, y; parallel=true)
        for i in eachindex(x, y)
            _x = x[i]
            _y = y[i]
            @test all(≈(f(_x, _y)), (itp(_x, _y), vals[i], vals2[i], vals3[i], vals4[i]))
        end
    end

    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30, add_ghost_triangles=true)
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(get_points(tri), z)
    xx = LinRange(0, 1, 50)
    yy = LinRange(0, 1, 50)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    for _ in 1:500
        vals = itp(x, y; parallel=false)
        vals2 = similar(vals)
        itp(vals2, x, y; parallel=false)
        vals3 = itp(x, y; parallel=true)
        vals4 = similar(vals3)
        itp(vals4, x, y; parallel=true)
        for i in eachindex(x, y)
            _x = x[i]
            _y = y[i]
            @test all(val -> isapprox(val, f(_x, _y), rtol=1e-1), (itp(_x, _y), vals[i], vals2[i], vals3[i], vals4[i]))
        end
    end
end


rng = StableRNG(123)
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 15, 15, add_ghost_triangles=true)
pts = get_points(tri)
z = [f(x, y) for (x, y) in each_point(tri)]
itp = interpolate(get_points(tri), z)
xx = LinRange(0, 1, 50)
yy = LinRange(0, 1, 50)
zz = [itp(x, y) for x in xx, y in yy]
_zz = [f(x, y) for x in xx, y in yy]

fig = Figure()
ax = Axis(fig[1, 1], aspect=1)
contourf!(ax, xx, yy, _zz, colormap=:viridis, levels=20)
ax = Axis(fig[1, 2], aspect=1)
contourf!(ax, xx, yy, zz, colormap=:viridis, levels=20)
fig