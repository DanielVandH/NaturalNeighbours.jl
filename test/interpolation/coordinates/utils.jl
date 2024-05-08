using ..NaturalNeighbours
using Test
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using Combinatorics
using StatsBase
using LinearAlgebra

@testset "Computing the Bowyer-Watson envelope" begin
    pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
    tri = triangulate(pts, randomise=false, delete_ghosts=false)
    envelope = Int64[]
    history = DT.InsertionEventHistory(tri)
    A = DT.Adjacent{Int64,NTuple{2,Int64}}()
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (6.0, 4.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 4.0)))
    @test DT.circular_equality(envelope, [1, 5, 7, 6, 1])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (6.0, 1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 1.0)))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (7.0, -1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (7.0, -1.0)))
    @test DT.circular_equality(envelope, [2, DT.GhostVertex, 3, 8, 7, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope(tri, (7.0, -1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (7.0, -1.0)))
    @test DT.circular_equality(envelope, [2, DT.GhostVertex, 3, 8, 7, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, (0.0, 3.0))
    @test DT.is_on(DT.point_position_relative_to_triangle(tri, V, (0.0, 3.0)))
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (6.0, 4.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 4.0)))
    for T in DT.each_added_triangle(history)
        i, j, k = triangle_vertices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A.adjacent) == 3length(DT.each_added_triangle(history))
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (6.0, 1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (6.0, 1.0)))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    for T in DT.each_added_triangle(history)
        i, j, k = triangle_vertices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A.adjacent) == 3length(DT.each_added_triangle(history))
    @test DT.circular_equality(envelope, [2, 3, 8, 7, 5, 2])
    envelope, A, history, V = NNI.compute_bowyer_envelope!(envelope, tri, history, A, (7.0, -1.0))
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, (7.0, -1.0)))
    for T in DT.each_added_triangle(history)
        i, j, k = triangle_vertices(T)
        @test DT.get_adjacent(A, i, j) == k
        @test DT.get_adjacent(A, j, k) == i
        @test DT.get_adjacent(A, k, i) == j
    end
    @test length(A.adjacent) == 3length(DT.each_added_triangle(history))
    @test DT.circular_equality(envelope, [2, DT.GhostVertex, 3, 8, 7, 2])
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
    boundary_nodes = get_convex_hull_vertices(tri)
    A = NNI.polygon_area(pts[:, boundary_nodes])
    _A = DT.get_area(tri)
    @test A ≈ _A
end

@testset "Testing sort_five" begin
    @test NNI.sort_five(5, 4, 3, 2, 1) == (1, 2, 3, 4, 5)
    @test NNI.sort_five(1, 2, 3, 4, 5) == (1, 2, 3, 4, 5)
    @test NNI.sort_five(3, 1, 4, 5, 2) == (1, 2, 3, 4, 5)
    @test NNI.sort_five(5, 5, 5, 5, 5) == (5, 5, 5, 5, 5)
    @test NNI.sort_five(1, 1, 2, 2, 2) == (1, 1, 2, 2, 2)
    @test NNI.sort_five(2, 2, 2, 1, 1) == (1, 1, 2, 2, 2)

    @testset "random tests" begin
        for _ in 1:10000
            arr = rand(1:5, 5)
            @test NNI.sort_five(arr...) == Tuple(sort(arr))
        end
    end
end

@testset "Testing count_unique_sorted" begin
    for _ in 1:10000
        arr = rand(1:5, 5)
        i, j, k, ℓ, m = NNI.sort_five(arr...)
        ct = NNI.count_unique_sorted(i, j, k, ℓ, m)
        @test ct == length(unique(arr))
    end
end

@testset "Testing group_sort" begin
    _up(i, j, k, ℓ, m) = (unique ∘ permutations)((i, j, k, ℓ, m))
    for _ in 1:10000
        arr = rand(1:5, 5)
        (a, b, c, d, e), case = NNI.group_sort(arr...)
        @inferred NNI.group_sort(arr...)
        if case == 1
            # iiiii
            @test a == b == c == d == e
        elseif case == 2
            # iiiij 
            @test (a == b == c == d) && (d ≠ e)
        elseif case == 3
            # iiijj 
            @test (a == b == c) && (c ≠ d) && (d == e)
        elseif case == 4
            # iiijk 
            @test (a == b == c) && (c ≠ d) && (d ≠ e)
        elseif case == 5
            # iijjk 
            @test (a == b) && (b ≠ c) && (c == d) && (d ≠ e) && (e ≠ a)
        elseif case == 6
            # iijkℓ 
            @test (a == b) && (b ≠ c) && (c ≠ d) && (d ≠ c) && (d ≠ a) && (e ≠ d) && (e ≠ c) && (e ≠ a)
        elseif case == 7
            # ijkℓm 
            @test (a ≠ b) && (b ≠ c) && (c ≠ d) && (d ≠ e) && (e ≠ a)
        end
    end
    for i in 1:5, perm in _up(i, i, i, i, i)
        _, case = NNI.group_sort(perm...)
        @test case == 1
    end
    for i in 1:5, j in 6:10, perm in _up(i, i, i, i, j)
        _, case = NNI.group_sort(perm...)
        @test case == 2
    end
    for i in 1:5, j in 6:10, perm in _up(i, i, i, j, j)
        _, case = NNI.group_sort(perm...)
        @test case == 3
    end
    for i in 1:5, j in 6:10, k in 11:15, perm in _up(i, i, i, j, k)
        _, case = NNI.group_sort(perm...)
        @test case == 4
    end
    for i in 1:5, j in 6:10, k in 11:15, perm in _up(i, i, j, j, k)
        _, case = NNI.group_sort(perm...)
        @test case == 5
    end
    for i in 1:5, j in 6:10, k in 11:15, ℓ in 16:20, perm in _up(i, i, j, k, ℓ)
        _, case = NNI.group_sort(perm...)
        @test case == 6
    end
    for i in 1:5, j in 6:10, k in 11:15, ℓ in 16:20, m in 21:25, perm in _up(i, j, k, ℓ, m)
        _, case = NNI.group_sort(perm...)
        @test case == 7
    end
end

@testset "directional_derivative" begin
    for _ in 1:200
        tri = triangulate(rand(2, 500))
        i = rand(1:10)
        j = rand(1:10)
        ∇ = [rand(2) for _ in 1:500]
        N₀ = sample(1:500, 10, replace=false)
        @test NNI.directional_derivative(tri, i, j, N₀, ∇) ≈ dot(get_point(tri, N₀[j]) .- get_point(tri, N₀[i]), ∇[N₀[i]])
    end
end

to_mat(H) = [H[1] H[3]; H[3] H[2]]
@testset "hessian_form" begin
    for _ in 1:200
        tri = triangulate(rand(2, 500))
        i = rand(1:10)
        j = rand(1:10)
        k = rand(1:10)
        H = [rand(3) for _ in 1:500]
        N₀ = sample(1:500, 10, replace=false)
        @test NNI.hessian_form(tri, i, j, k, N₀, H) ≈ collect(get_point(tri, N₀[j]) .- get_point(tri, N₀[i]))' * to_mat(H[N₀[i]]) * collect(get_point(tri, N₀[k]) .- get_point(tri, N₀[i]))
    end
end