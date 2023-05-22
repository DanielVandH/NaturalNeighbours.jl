using ..NaturalNeighbours
using Test
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation

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
    envelope, A, history, V = NNI.compute_bowyer_envelope(tri, (7.0, -1.0))
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