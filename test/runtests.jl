using NaturalNeighbours
using Test

@testset "NaturalNeighbours.jl" begin
    # Write your tests here.
end

using DelaunayTriangulation, Random, CairoMakie, StableRNGs, LinearAlgebra, ReferenceTests, ElasticArrays, Optimization, OptimizationNLopt, StatsBase
const DT = DelaunayTriangulation
const NNI = NaturalNeighbours

include("helper_functions/point_generator.jl")
include("helper_functions/slow_derivative_tests.jl")
include("helper_functions/test_functions.jl")

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

@testset "Natural coordinates" begin
    for method in (:sibson, :triangle, :nearest, :laplace)
        method = NNI.iwrap(method)
        pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
        tri = triangulate(pts, randomise=false, delete_ghosts=false)
        n = 2500
        pts = random_points_in_convex_hull(tri, n)
        for p in Iterators.flatten((pts, each_point(tri)))
            natural_coordinates = NNI.compute_natural_coordinates(method, tri, p)
            @test sum(NNI.get_coordinates(natural_coordinates)) ≈ 1
            δ = NNI.get_barycentric_deviation(natural_coordinates)
            if method ≠ NNI.Nearest()
                @test δ ≈ 0 atol = 1e-5
            else
                @test δ ≈ norm(p .- get_point(tri, NNI.get_indices(natural_coordinates)[1]))
            end
        end

        for _ in 1:50
            pts = [(randn() + rand(), rand() - 2randn()) for _ in 1:250]
            tri = triangulate(pts, delete_ghosts=false)
            n = 5000
            random_points = random_points_in_convex_hull(tri, n)
            for p in Iterators.flatten((random_points, each_point(tri)))
                natural_coordinates = NNI.compute_natural_coordinates(method, tri, p)
                @test sum(NNI.get_coordinates(natural_coordinates)) ≈ 1
                δ = NNI.get_barycentric_deviation(natural_coordinates)
                if method ≠ NNI.Nearest()
                    @test δ ≈ 0 atol = 1e-5
                else
                    @test δ ≈ norm(p .- get_point(tri, NNI.get_indices(natural_coordinates)[1]))
                end
            end
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

@testset "iterated_neighbourhood!" begin
    rng = StableRNG(123)
    for _ in 1:50
        pts = [(rand(rng), rand(rng)) for _ in 1:50]
        tri = triangulate(pts, rng=rng, delete_ghosts=false)
        f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
        z = [f(x, y) for (x, y) in pts]
        S = Set{Int64}()
        S′ = Set{Int64}()
        n_cache = NNI.NaturalNeighboursCache(tri)
        i = 7
        d = 1
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, i, d, n_cache)
        @test all(isone, λ)
        @test sort(collect(E)) == sort(collect(DT.iterated_neighbourhood(tri, i, 1)))
        p = get_point(tri, i)
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache)
        @test all(isone, λ)
        @test λ == 1
        @test sort(collect(E)) == sort(collect(DT.iterated_neighbourhood(tri, i, 1)))
        p = random_points_in_convex_hull(tri, 1)[1]
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, NNI.NaturalNeighboursCache(tri); rng=StableRNG(881))
        nc = NNI.compute_natural_coordinates(Sibson(), tri, p, NNI.NaturalNeighboursCache(tri); rng=StableRNG(881))
        @test λ == nc.coordinates
        @test E == nc.indices
        p = random_points_in_convex_hull(tri, 1)[1]
        for _ in 1:100 # test rng is passed correctly
            λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, NNI.NaturalNeighboursCache(tri); rng=StableRNG(125))
            nc = NNI.compute_natural_coordinates(Sibson(), tri, p, NNI.NaturalNeighboursCache(tri); rng=StableRNG(125))
            @test λ == nc.coordinates
            @test E == nc.indices
        end

        d = 2
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, i, d, n_cache)
        @test λ == 1
        _S = DT.iterated_neighbourhood(tri, i, 2)
        @test sort(collect(E)) == sort(collect(_S))
        p = get_point(tri, i)
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache)
        @test λ == 1
        @test sort(collect(E)) == sort(collect(_S))
        p = random_points_in_convex_hull(tri, 1)[1]
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache)
        nc = NNI.compute_natural_coordinates(Sibson(), tri, p)
        _S = [get_neighbours(tri, i) for i in nc.indices]
        _S1 = copy(nc.indices)
        push!(_S1, reduce(union, _S)...)
        filter!(!DT.is_boundary_index, _S1)
        unique!(_S1)
        @test sort(E) == sort(_S1)
        for i in eachindex(E)
            if 1 ≤ i ≤ length(λ)
                @test NNI.get_λ(λ, i, true) == λ[i]
                @test NNI.get_λ(λ, i, false) == 1
            else
                @test NNI.get_λ(λ, i, true) == NNI.get_λ(λ, i, false) == 1
            end
        end

        d = 3
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, i, d, n_cache)
        @test λ == 1
        _S = DT.iterated_neighbourhood(tri, i, 3)
        @test sort(collect(E)) == sort(collect(_S))
        p = get_point(tri, i)
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache)
        @test λ == 1
        @test sort(collect(E)) == sort(collect(_S))
        p = random_points_in_convex_hull(tri, 1)[1]
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache)
        nc = NNI.compute_natural_coordinates(Sibson(), tri, p)
        _S = [get_neighbours(tri, i) for i in nc.indices]
        _S1 = copy(nc.indices)
        push!(_S1, reduce(union, _S)...)
        filter!(!DT.is_boundary_index, _S1)
        unique!(_S1)
        _S2 = [get_neighbours(tri, i) for i in _S1]
        push!(_S1, reduce(union, _S2)...)
        filter!(!DT.is_boundary_index, _S1)
        unique!(_S1)
        @test sort(E) == sort(_S1)
        for i in eachindex(E)
            if 1 ≤ i ≤ length(λ)
                @test NNI.get_λ(λ, i, true) == λ[i]
                @test NNI.get_λ(λ, i, false) == 1
            else
                @test NNI.get_λ(λ, i, true) == NNI.get_λ(λ, i, false) == 1
            end
        end
    end
end

@testset "Interpolation" begin
    rng = StableRNG(123)
    pts = [(rand(rng), rand(rng)) for _ in 1:50]
    tri = triangulate(pts, rng=rng, delete_ghosts=false)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    z = [f(x, y) for (x, y) in pts]

    itp = interpolate(tri, z; derivatives=true, parallel=false)
    @test DT.get_triangulation(itp) == tri
    @test NNI.get_z(itp) == z
    @test length(NNI.get_neighbour_cache(itp)) == Base.Threads.nthreads()
    @test length(NNI.get_derivative_cache(itp)) == Base.Threads.nthreads()
    @test NNI.get_neighbour_cache(itp, 1) == itp.neighbour_cache[1]
    @test NNI.get_neighbour_cache(itp, 2) == itp.neighbour_cache[2]
    @test NNI.get_derivative_cache(itp) == itp.derivative_cache
    @test NNI.get_derivative_cache(itp, 1) == itp.derivative_cache[1]
    @test NNI.get_derivative_cache(itp, 2) == itp.derivative_cache[2]
    @test NNI.get_gradient(itp) == itp.gradient
    @test !isnothing(NNI.get_gradient(itp))
    @test NNI.get_gradient(itp, 1) == itp.gradient[1]
    @test NNI.get_gradient(itp, 2) == itp.gradient[2]
    @test NNI.get_hessian(itp) == itp.hessian
    @test !isnothing(NNI.get_hessian(itp))
    @test NNI.get_hessian(itp, 1) == itp.hessian[1]
    @test NNI.get_hessian(itp, 2) == itp.hessian[2]
    _itp =  interpolate(tri, z; derivatives=false, parallel=false)
    @test NNI.get_gradient(_itp) === nothing
    @test NNI.get_hessian(_itp) === nothing
    @test itp isa NNI.NaturalNeighboursInterpolant
    DT.lock_convex_hull!(tri)
    @test_throws ArgumentError interpolate(tri, z)
    DT.unlock_convex_hull!(tri)
    @test_throws AssertionError interpolate(tri, z[1:end-1])
    w =rand(length(z))
    y = rand(length(z))
    __itp =  interpolate(tri, z; derivatives=false, parallel=false, gradient = w)
    @test NNI.get_gradient(__itp) === w
    __itp =  interpolate(tri, z; derivatives=false, parallel=false, gradient = w, hessian = y)
    @test NNI.get_gradient(__itp) === w
    @test NNI.get_hessian(__itp) === y

    x = getx.(pts)
    y = gety.(pts)
    test_interpolant(itp, x, y, f)
    test_interpolant(itp, x, y, z)

    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30, add_ghost_triangles=true)
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(get_points(tri), z; derivatives=true)
    xx = LinRange(0, 1, 50)
    yy = LinRange(0, 1, 50)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    test_interpolant(itp, x, y, f)
    x = getx.(get_points(tri))
    y = gety.(get_points(tri))
    test_interpolant(itp, x, y, f)
    test_interpolant(itp, x, y, z)
end

@testset "Basic extrapolation" begin
    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 5, 10, add_ghost_triangles=true)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    pts = get_points(tri)
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(pts, z)

    p = (1.5, 0.7)
    V = jump_and_march(tri, p)
    _V = DT.rotate_ghost_triangle_to_standard_form(V)
    i, j = indices(_V)
    a, b = get_point(tri, i, j)
    dab = norm(b .- a)
    dbp = norm((1.0, 0.7) .- b)
    t = dbp / dab
    _z = t * z[i] + (1 - t) * z[j]
    __z = itp(getx(p), gety(p); method=:triangle)
    @test _z ≈ itp(getx(p), gety(p); method=:triangle)
    @test _z ≈ itp(getx(p), gety(p); method=:sibson)
    @test _z ≈ itp(getx(p), gety(p); method=:laplace)
    @test __z ≈ itp(1.8, 0.7; method=:triangle)
    @test __z ≈ itp(1.8, 0.7; method=:sibson)
    @test __z ≈ itp(1.8, 0.7; method=:laplace)
end

@testset "Example I: No extrapolation" begin
    ## Example I: No extrapolation 
    # Define the interpolant
    rng = StableRNG(123)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    x = vec([(i - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
    y = vec([(j - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
    z = f.(x, y)
    itp = interpolate(x, y, z; rng)

    # Points to evaluate the interpolant at 
    xx = LinRange(0, 1, 50)
    yy = LinRange(0, 1, 50)
    _x = vec([x for x in xx, _ in yy])
    _y = vec([y for _ in xx, y in yy])

    # Evaluate the interpolant
    sibson_vals = itp(_x, _y; method=:sibson, rng) # multithreaded
    triangle_vals = itp(_x, _y; method=:triangle, rng)
    exact_vals = [f(x, y) for x in xx, y in yy]
    sibson_vals = reshape(sibson_vals, (length(xx), length(yy)))
    triangle_vals = reshape(triangle_vals, (length(xx), length(yy)))

    # Get the errors 
    sibson_errs = abs.(sibson_vals .- exact_vals) ./ abs.(exact_vals)
    triangle_errs = abs.(triangle_vals .- exact_vals) ./ abs.(exact_vals)

    # Plot the results
    fig = Figure(fontsize=33)
    make_ax = (i, j, title) -> begin
        Axis(fig[i, j], title=title, titlealign=:left,
            width=400, height=400,
            xticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]), yticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]))
    end
    ax1 = make_ax(1, 1, L"(a):$ $ Sibson")
    contourf!(ax1, xx, yy, sibson_vals, colormap=:viridis, levels=-1:0.05:0)
    ax2 = make_ax(1, 2, L"(b):$ $ Triangle")
    contourf!(ax2, xx, yy, triangle_vals, colormap=:viridis, levels=-1:0.05:0)
    ax3 = make_ax(1, 3, L"(c):$ $ Exact")
    contourf!(ax3, xx, yy, exact_vals, colormap=:viridis, levels=-1:0.05:0)
    ax4 = make_ax(2, 3, L"(f):$ $ Data")
    tricontourf!(ax4, x, y, z, colormap=:viridis, levels=-1:0.05:0)
    ax5 = make_ax(2, 1, L"(d):$ $ Sibson error")
    contourf!(ax5, xx, yy, sibson_errs, colormap=:viridis, levels=0:0.01:0.1)
    ax6 = make_ax(2, 2, L"(e):$ $ Triangle error")
    contourf!(ax6, xx, yy, triangle_errs, colormap=:viridis, levels=0:0.01:0.1)
    for ax in (ax1, ax2, ax3, ax3, ax4, ax5, ax6)
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
        scatter!(ax, x, y, markersize=9, color=:red)
    end
    resize_to_layout!(fig)
    @test_reference "figures/example_1.png" fig
end

@testset "Example II: Extrapolation" begin
    ## Define the interpolant 
    rng = StableRNG(1235)
    f = (x, y) -> (1 / 9) * (tanh(9y - 9x) + 1)
    x = rand(rng, 100)
    y = rand(rng, 100)
    z = f.(x, y)
    itp = interpolate(x, y, z; rng)

    ## Points to evaluate at 
    xx = LinRange(-1 / 2, 3 / 2, 250)
    yy = LinRange(-1 / 2, 3 / 2, 250)
    _x = vec([x for x in xx, _ in yy])
    _y = vec([y for _ in xx, y in yy])

    ## Evaluate the interpolant
    sibson_vals = itp(_x, _y; method=:sibson) # multithreaded
    triangle_vals = itp(_x, _y; method=:triangle)
    exact_vals = [f(x, y) for x in xx, y in yy]
    sibson_vals = reshape(sibson_vals, (length(xx), length(yy)))
    triangle_vals = reshape(triangle_vals, (length(xx), length(yy)))
    sibson_errs = abs.(sibson_vals .- exact_vals) ./ (1.0 .+ abs.(exact_vals))
    triangle_errs = abs.(triangle_vals .- exact_vals) ./ (1.0 .+ abs.(exact_vals))

    ## Plot 
    fig = Figure(fontsize=33)
    make_ax = (i, j, title) -> begin
        Axis(fig[i, j], title=title, titlealign=:left,
            width=400, height=400,
            xticks=([-0.5, 0, 0.5, 1.0, 1.5], [L"-0.5", L"0", L"0.5", L"1", L"1.5"]),
            yticks=([-0.5, 0, 0.5, 1.0, 1.5], [L"-0.5", L"0", L"0.5", L"1", L"1.5"]))
    end
    ax1 = make_ax(1, 1, L"(a):$ $ Sibson")
    contourf!(ax1, xx, yy, sibson_vals, colormap=:viridis, levels=-0.4:0.1:0.4)
    ax2 = make_ax(1, 2, L"(b):$ $ Triangle")
    contourf!(ax2, xx, yy, triangle_vals, colormap=:viridis, levels=-0.4:0.1:0.4)
    ax3 = make_ax(1, 3, L"(c):$ $ Exact")
    contourf!(ax3, xx, yy, exact_vals, colormap=:viridis, levels=-0.4:0.1:0.4)
    ax4 = make_ax(2, 3, L"(f):$ $ Data")
    tricontourf!(ax4, x, y, z, colormap=:viridis, levels=-0.4:0.1:0.4)
    ax5 = make_ax(2, 1, L"(d):$ $ Sibson error")
    contourf!(ax5, xx, yy, sibson_errs, colormap=:viridis, levels=0:0.025:0.5)
    ax6 = make_ax(2, 2, L"(e):$ $ Triangle error")
    contourf!(ax6, xx, yy, triangle_errs, colormap=:viridis, levels=0:0.025:0.5)
    tri = itp.triangulation
    ch_idx = get_convex_hull_indices(tri)
    ch_points = [get_point(tri, i) for i in ch_idx]
    for ax in (ax1, ax2, ax3, ax3, ax4, ax5, ax6)
        xlims!(ax, -1 / 2, 3 / 2)
        ylims!(ax, -1 / 2, 3 / 2)
        lines!(ax, ch_points, color=:black, linewidth=3)
        scatter!(ax, x, y, markersize=9, color=:red)
    end
    resize_to_layout!(fig)
    @test_reference "figures/example_2.png" fig
end

@testset "Test coefficient values for each method" begin # used GeoGebra
    for _ in 1:100 # make sure rng is getting passed consistently 
        # Setup
        rng = StableRNG(872973)
        pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0),
            (14.0, 8.0), (4.0, 4.0), (10.0, 6.0),
            (6.0, 2.0), (12.0, 4.0), (0.0, 4.0),
            (2.5, 5.0), (7.0, 3.3), (4.5, 5.2),
            (13.0, 0.5), (12.0, 6.0), (8.5, 3.5),
            (0.5, 6.0), (1.5, 6.0), (3.5, 6.0),
            (0.5, 2.0), (2.5, 2.0), (2.5, 2.5),
            (9.0, 2.0), (8.5, 6.0), (4.0, 2.0)]
        tri = triangulate(pts, randomise=false, delete_ghosts=false, rng=rng)
        vorn = voronoi(tri, false)
        q = (5.0, 4.0)
        tri2 = deepcopy(tri)
        add_point!(tri2, q, rng=rng)
        vorn2 = voronoi(tri2, false)
        V = get_polygon(vorn2, num_points(tri2))
        AX2 = get_area(vorn2, num_points(tri2))

        # Sibson
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q; rng=rng)
        AF1G1D1B1A1 = 1.1739978952813 # E = 5
        A1B1C1 = 0.062500000375 # Z = 24
        B1D1E1C1 = 0.2540749084958 # G = 7
        H1I1E1D1G1 = 0.7376579777378 # K = 11
        F1G1J1K1 = 0.9489911839831 # L = 12
        K1I1J1 = 0.003313286868 # W = 23
        AX = AF1G1D1B1A1 + A1B1C1 + B1D1E1C1 + H1I1E1D1G1 + F1G1J1K1 + K1I1J1
        @test AX ≈ AX2 rtol = 1e-3
        @test nc.indices == [23, 12, 5, 24, 7, 11]
        @test nc.coordinates ≈ [K1I1J1, F1G1J1K1, AF1G1D1B1A1, A1B1C1, B1D1E1C1, H1I1E1D1G1] ./ AX rtol = 1e-2

        # Triangle 
        nc = NNI.compute_natural_coordinates(NNI.Triangle(), tri, q; rng=rng)
        V = jump_and_march(tri, q; rng)
        @test nc.indices == [5, 11, 12]
        @test nc.coordinates ≈ [0.52, 0.3, 0.18] rtol = 1e-2

        # Laplace 
        nc = NNI.compute_natural_coordinates(NNI.Laplace(), tri, q; rng=rng)
        dqw = 4.0311288741493
        dqk = 2.1189620100417
        dqg = 2.2360679774998
        dqz = 2.2360679774998
        dqe = 1.0
        dqℓ = 1.3
        dc1b1 = 2.2893491697301
        da1b1 = 0.0572608008105
        df1b1 = 2.2260773066834
        df1e1 = 1.4888476232694
        de1d1 = 0.5650898843856
        dd1c1 = 0.9335156761474
        k = dc1b1 / dqk
        w = da1b1 / dqw
        ℓ = df1b1 / dqℓ
        e = df1e1 / dqe
        z = de1d1 / dqz
        g = dd1c1 / dqg
        tot = k + w + ℓ + e + z + g
        k /= tot
        w /= tot
        ℓ /= tot
        e /= tot
        z /= tot
        g /= tot
        @test nc.indices == [23, 12, 5, 24, 7, 11]
        @test nc.coordinates ≈ [w, ℓ, e, z, g, k] rtol = 1e-2

        # Nearest 
        nc = NNI.compute_natural_coordinates(NNI.Nearest(), tri, q; rng=rng)
        @test nc.indices == [5]
        @test nc.coordinates ≈ [1.0]

        # Sibson(1)
        tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
        tri = triangulate(get_points(tri), randomise=false)
        f = (x, y) -> sin(x - y) + cos(x + y)
        z = [f(x, y) for (x, y) in each_point(tri)]
        itp = interpolate(tri, z; derivatives=true)
        q = (5.0, 5.0) # a data site
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q; rng=rng)
        ∇ = NNI.get_gradient(itp)
        ζ, α, β = NNI._compute_sibson_1_coordinates(nc, tri, z, ∇)
        @test ζ == 0.0
        @test α == 1.0
        @test β == 0.0
        q = (5.37841, 1.3881)
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q; rng=rng)
        ζ1, α1, β1 = NNI._compute_sibson_1_coordinates(nc, tri, z, ∇)
        λ, N₀ = NNI.get_coordinates(nc), NNI.get_indices(nc)
        r = [norm(q .- get_point(tri, i)) for i in N₀]
        γ = λ ./ r
        ζ = z[N₀] .+ [dot(∇[i], q .- get_point(tri, i)) for i in N₀]
        α = dot(λ, r) / sum(γ)
        β = dot(λ, r .^ 2)
        ζ = dot(ζ, γ) / sum(γ)
        @test α ≈ α1
        @test β ≈ β1
        @test ζ ≈ ζ1
    end
end

@testset "DerivativeCache" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    derivative_cache = NNI.DerivativeCache(tri)
    @test NNI.get_iterated_neighbourhood(derivative_cache) == derivative_cache.iterated_neighbourhood == Set{Int64}()
    @test NNI.get_second_iterated_neighbourhood(derivative_cache) == derivative_cache.second_iterated_neighbourhood == Set{Int64}()
    @test NNI.get_linear_matrix(derivative_cache) == derivative_cache.linear_matrix == ElasticMatrix{Float64}(undef, 2, 0)
    @test NNI.get_quadratic_matrix(derivative_cache) == derivative_cache.quadratic_matrix == ElasticMatrix{Float64}(undef, 9, 0)
    @test NNI.get_rhs_vector(derivative_cache) == derivative_cache.rhs_vector == Float64[]
    @test NNI.get_linear_sol(derivative_cache) == derivative_cache.linear_sol == [0.0, 0.0]
    @test NNI.get_quadratic_sol(derivative_cache) == derivative_cache.quadratic_sol == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test NNI.get_quadratic_matrix_no_cubic(derivative_cache) == derivative_cache.quadratic_matrix_no_cubic == ElasticMatrix{Float64}(undef, 5, 0)
    @test NNI.get_quadratic_sol_no_cubic(derivative_cache) == derivative_cache.quadratic_sol_no_cubic == [0.0, 0.0, 0.0, 0.0, 0.0]
end

@testset "iwrap" begin
    @test NNI.iwrap(NNI.Sibson()) == NNI.Sibson()
    @test NNI.iwrap(NNI.Triangle()) == NNI.Triangle()
    @test NNI.iwrap(NNI.Nearest()) == NNI.Nearest()
    @test NNI.iwrap(NNI.Laplace()) == NNI.Laplace()
    @test NNI.iwrap(:sibson) == NNI.Sibson()
    @test NNI.iwrap(:triangle) == NNI.Triangle()
    @test NNI.iwrap(:nearest) == NNI.Nearest()
    @test NNI.iwrap(:laplace) == NNI.Laplace()
    @test_throws ArgumentError NNI.iwrap(:lap)

    @test NNI.iwrap(NNI.Sibson(1)) == NNI.Sibson(1)
    @test NNI.iwrap(NNI.Triangle(1)) == NNI.Triangle(0)
    @test_throws ArgumentError NNI.Sibson(5)
    @test NNI.iwrap(NNI.Laplace(1)) == NNI.Laplace(0)
end

@testset "dwrap" begin
    @test NNI.dwrap(NNI.Direct()) == NNI.Direct()
    @test NNI.dwrap(:direct) == NNI.Direct()
    @test_throws ArgumentError NNI.dwrap(:dir)
    @test NNI.dwrap(NNI.Iterative()) == NNI.Iterative()
    @test NNI.dwrap(:iterative) == NNI.Iterative()
    @test_throws ArgumentError NNI.dwrap(:iter)
end

@testset "Estimating derivatives with weighted least squares" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> sin(x - y) + cos(x + y)
    f′ = (x, y) -> [cos(x - y) - sin(x + y), -cos(x - y) - sin(x + y)]
    f′′ = (x, y) -> [-sin(x - y)-cos(x + y) sin(x - y)-cos(x + y)
        sin(x - y)-cos(x + y) -sin(x - y)-cos(x + y)]
    z = [f(x, y) for (x, y) in each_point(tri)]
    @testset "Direct" begin
        @testset "At data sites" begin
            flag = 0
            for _ in 1:100
                i = rand(1:num_points(tri))
                p = get_point(tri, i)

                # Gradient
                ∇opt, ∇ls = estimate_gradient_direct(tri, i, z)
                @test ∇opt ≈ ∇ls rtol = 1e-2
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, i, 1)
                ∇2 = (
                    NNI.generate_first_order_derivatives(NNI.Direct(), tri, z, z[i], i, λ, E),
                    NNI._generate_first_order_derivatives_direct(tri, z, z[i], i, λ, E)
                )
                @test collect(∇2[1]) ≈ collect(∇2[2]) # not exactly == sometimes due to rng
                @test collect(∇2[1]) ≈ ∇opt rtol = 1e-3
                flag += isapprox(f′(p...), collect(∇2[1]), rtol=1e-1)

                # Hessian: Cubic
                (∇opt, ℋopt), (∇ls, ℋls) = estimate_gradient_hessian_cubic_direct(tri, i, z)
                @test ∇opt ≈ ∇ls rtol = 1e-4
                @test ℋopt ≈ ℋls rtol = 1e-4
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, i, 3)
                ∇ℋ2 = (
                    NNI.generate_second_order_derivatives(NNI.Direct(), tri, z, z[i], i, λ, E),
                    NNI._generate_second_order_derivatives_direct(tri, z, z[i], i, E, NNI.DerivativeCache(tri))
                )
                ∇2_1, ℋ2_1 = ∇ℋ2[1]
                ∇2_2, ℋ2_2 = ∇ℋ2[2]
                @test collect(∇2_1) ≈ collect(∇2_2)
                @test collect(ℋ2_1) ≈ collect(ℋ2_2)
                @test collect(∇2_1) ≈ ∇opt rtol = 1e-4
                @test collect(ℋ2_1) ≈ ℋopt rtol = 1e-4
                @test f′(p...) ≈ collect(∇2_1) rtol = 1e-2
                @test f′′(p...) ≈ [ℋ2_1[1] ℋ2_1[3]; ℋ2_1[3] ℋ2_1[2]] rtol = 1e-1

                # Hessian: Quadratic 
                (∇opt, ℋopt), (∇ls, ℋls) = estimate_gradient_hessian_quadratic_direct(tri, i, z)
                @test ∇opt ≈ ∇ls rtol = 1e-4
                @test ℋopt ≈ ℋls rtol = 1e-4
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, i, 2)
                ∇ℋ2 = (
                    NNI.generate_second_order_derivatives(NNI.Direct(), tri, z, z[i], i, λ, E; use_cubic_terms=false),
                    NNI._generate_second_order_derivatives_direct(tri, z, z[i], i, E; use_cubic_terms=false)
                )
                ∇2_1, ℋ2_1 = ∇ℋ2[1]
                ∇2_2, ℋ2_2 = ∇ℋ2[2]
                @test collect(∇2_1) ≈ collect(∇2_2)
                @test collect(ℋ2_1) ≈ collect(ℋ2_2)
                @test collect(∇2_1) ≈ ∇opt rtol = 1e-5
                @test collect(ℋ2_1) ≈ ℋopt rtol = 1e-5
                @test f′(p...) ≈ collect(∇2_1) rtol = 1e-1
                @test f′′(p...) ≈ [ℋ2_1[1] ℋ2_1[3]; ℋ2_1[3] ℋ2_1[2]] rtol = 0.5 atol = 1e-1
            end
            @test flag / 100 > 0.95
        end

        @testset "At off-site points" begin
            flag = 0
            rng = StableRNG(35)
            for _ in 1:100
                # Gradient
                itp = interpolate(tri, z; derivatives=false)
                p = random_points_in_convex_hull(tri, 1; rng)[1]
                ∇opt, ∇ls = estimate_gradient_direct(tri, p, z)
                @test ∇opt ≈ ∇ls rtol = 1e-1
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, p, 1)
                ∇2 = (
                    NNI.generate_first_order_derivatives(NNI.Direct(), tri, z, itp(p...), p, λ, E),
                    NNI._generate_first_order_derivatives_direct(tri, z, itp(p...), p, λ, E)
                )
                @test collect(∇2[1]) ≈ collect(∇2[2]) rtol = 1e-1
                @test collect(∇2[1]) ≈ ∇opt rtol = 0.3
                @test f′(p...) ≈ collect(∇2[1]) rtol = 1e-1 atol = 0.2

                # Hessian: Cubic 
                p = random_points_in_convex_hull(tri, 1; rng)[1]
                (∇opt, ℋopt), (∇ls, ℋls) = estimate_gradient_hessian_cubic_direct(tri, p, z)
                @test ∇opt ≈ ∇ls rtol = 1e-2
                @test ℋopt ≈ ℋls rtol = 1e-2
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, p, 3)
                ∇ℋ2 = (
                    NNI.generate_second_order_derivatives(NNI.Direct(), tri, z, itp(p...), p, λ, E),
                    NNI._generate_second_order_derivatives_direct(tri, z, itp(p...), p, E, NNI.DerivativeCache(tri))
                )
                ∇2_1, ℋ2_1 = ∇ℋ2[1]
                ∇2_2, ℋ2_2 = ∇ℋ2[2]
                @test collect(∇2_1) ≈ collect(∇2_2)
                @test collect(ℋ2_1) ≈ collect(ℋ2_2)
                @test collect(∇2_1) ≈ ∇opt rtol = 1e-4
                @test collect(ℋ2_1) ≈ ℋopt rtol = 1e-4
                @test f′(p...) ≈ collect(∇2_1) rtol = 1e-1
                flag += isapprox(f′′(p...), [ℋ2_1[1] ℋ2_1[3]; ℋ2_1[3] ℋ2_1[2]], rtol=1e-1, atol=1e-1)

                # Hessian: Quadratic 
                (∇opt, ℋopt), (∇ls, ℋls) = estimate_gradient_hessian_quadratic_direct(tri, p, z)
                @test ∇opt ≈ ∇ls rtol = 1e-4
                @test ℋopt ≈ ℋls rtol = 1e-4
                λ, E = NNI.get_taylor_neighbourhood!(Set{Int64}(), Set{Int64}(), tri, p, 2)
                ∇ℋ2 = (
                    NNI.generate_second_order_derivatives(NNI.Direct(), tri, z, itp(p...), p, λ, E; use_cubic_terms=false),
                    NNI._generate_second_order_derivatives_direct(tri, z, itp(p...), p, E, NNI.DerivativeCache(tri); use_cubic_terms=false)
                )
                ∇2_1, ℋ2_1 = ∇ℋ2[1]
                ∇2_2, ℋ2_2 = ∇ℋ2[2]
                @test collect(∇2_1) ≈ collect(∇2_2)
                @test collect(ℋ2_1) ≈ collect(ℋ2_2)
                @test collect(∇2_1) ≈ ∇opt rtol = 1e-4
                @test collect(ℋ2_1) ≈ ℋopt rtol = 1e-4
                @test f′(p...) ≈ collect(∇2_1) rtol = 1e-1
                flag += isapprox(f′′(p...), [ℋ2_1[1] ℋ2_1[3]; ℋ2_1[3] ℋ2_1[2]], rtol=0.2)
            end
            @test flag / 200 > 0.8
        end
    end

    @testset "Iterative" begin
        tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
        tri = triangulate(get_points(tri), randomise=false)
        f = (x, y) -> sin(x - y) + cos(x + y)
        f′ = (x, y) -> [cos(x - y) - sin(x + y), -cos(x - y) - sin(x + y)]
        f′′ = (x, y) -> [-sin(x - y)-cos(x + y) sin(x - y)-cos(x + y)
            sin(x - y)-cos(x + y) -sin(x - y)-cos(x + y)]
        z = [f(x, y) for (x, y) in each_point(tri)]

        nt = Base.Threads.nthreads()
        derivative_caches = [NNI.DerivativeCache(tri) for _ in 1:nt]
        neighbour_caches = [NNI.NaturalNeighboursCache(tri) for _ in 1:nt]

        derivative_method = :iterative
        method = NNI.dwrap(derivative_method)
        parallel_derivatives = true
        initial_gradients = NNI.generate_gradients(tri, z, derivative_caches, neighbour_caches; method, parallel=parallel_derivatives)
        _initial_gradients = deepcopy(initial_gradients)
        for i in eachindex(initial_gradients)
            G1, G2 = estimate_gradient_direct(tri, i, z; use_sibson_weight=true)
            G3 = collect(initial_gradients[i])
            @test G1 ≈ G2 rtol = 1e-4
            @test G1 ≈ G3 rtol = 1e-4
            @test G2 ≈ G3 rtol = 1e-4
        end
        Gpar = NNI.generate_gradients(tri, z, derivative_caches, neighbour_caches; method, parallel=true)
        Gser = NNI.generate_gradients(tri, z, derivative_caches, neighbour_caches; method, parallel=false)
        @test Gpar == Gser
        ∇par, ℋpar = NNI.generate_derivatives(tri, z, derivative_caches, neighbour_caches; method, initial_gradients, parallel=true)
        ∇ser, ℋser = NNI.generate_derivatives(tri, z, derivative_caches, neighbour_caches; method, initial_gradients, parallel=false)
        @test initial_gradients == _initial_gradients # make sure initial_gradients is not modified
        @test ∇par == ∇ser
        @test ℋpar == ℋser
        _α = (0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5)
        flags = zeros(Int64, 8, length(_α))
        for (j, α) in enumerate(_α)
            ∇, ℋ = NNI.generate_derivatives(tri, z, derivative_caches, neighbour_caches; method, initial_gradients, parallel=true, alpha=α)
            for i in eachindex(initial_gradients)
                (G1, H1), (G2, H2) = estimate_gradient_hessian_from_initial_gradients(tri, i, z, α; initial_gradients)
                G3 = collect(∇[i])
                H3 = collect(ℋ[i])
                flags[1, j] += isapprox(G1, G2, rtol=1e-2)
                flags[2, j] += isapprox(G1, G3, rtol=1e-2)
                flags[3, j] += isapprox(G2, G3, rtol=1e-2)
                flags[4, j] += isapprox(H1, H2, rtol=1e-1)
                flags[5, j] += isapprox(H1, H3, rtol=1e-1)
                flags[6, j] += isapprox(H2, H3, rtol=1e-1)
                pᵢ = get_point(tri, i)
                xᵢ, yᵢ = getxy(pᵢ)
                G4 = f′(xᵢ, yᵢ)
                H4 = f′′(xᵢ, yᵢ)[[1, 4, 2]]
                flags[7, j] += isapprox(G3, G4, rtol=1e-1, atol=1e-1)
                flags[8, j] += isapprox(H3, H4, rtol=1e-1, atol=1e-1)
            end
        end
        normalised_flags = flags ./ length(initial_gradients)
        @test all(>(0.85), normalised_flags[:, 1])
        @test all(>(0.85), normalised_flags[:, 2])
        @test all(>(0.85), normalised_flags[:, 3])
        @test all(>(0.85), normalised_flags[:, 4])
        @test all(>(0.85), normalised_flags[:, 5])
        @test all(>(0.85), normalised_flags[:, 6])
        @test all(>(0.75), normalised_flags[:, 7])
        @test all(>(0.75), normalised_flags[[2, 4, 5, 6, 7, 8], 8])
        @test all(>(0.15), normalised_flags[[1, 3], 8])
    end
end

@testset "get_λ" begin
    @test NNI.get_λ([1.0, 2.0, 3.0], 2, true) == 2.0
    @test NNI.get_λ([1.0, 2.0, 3.0, 4.0], 5, true) == 1.0
    @test NNI.get_λ([2.3, 5.0], 1, false) == 1.0
    @test NNI.get_λ(1.0, 3, true) == 1.0
end

@testset "Checking that keyword arguments are passed fine" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(tri, z)
    ∂ = differentiate(itp, 2)
    @test_throws ArgumentError ∂(0.0, 0.0; method=Iterative())
    @test_throws ArgumentError ∂(0.0, 0.0; method=:iterative)
    itp1 = interpolate(tri, z; derivatives=false)
    itp2 = interpolate(tri, z; derivatives=true)
    ∂11 = differentiate(itp1, 1)
    ∂12 = differentiate(itp1, 2)
    ∂21 = differentiate(itp2, 1)
    ∂22 = differentiate(itp2, 2)

    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Direct(), interpolant_method=Sibson(),
        alpha=0.01, use_cubic_terms=true,
        use_sibson_weight=true,
        tri=tri, z=z)
    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Direct(), interpolant_method=Laplace(),
        alpha=0.01, use_cubic_terms=true,
        use_sibson_weight=true,
        tri=tri, z=z)
    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Direct(), interpolant_method=Laplace(),
        alpha=0.01, use_cubic_terms=true,
        use_sibson_weight=false,
        tri=tri, z=z)
    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Direct(), interpolant_method=Sibson(),
        alpha=0.01, use_cubic_terms=false,
        use_sibson_weight=true,
        tri=tri, z=z)
    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Direct(), interpolant_method=Sibson(),
        alpha=0.13, use_cubic_terms=false,
        use_sibson_weight=true,
        tri=tri, z=z)
    @test_throws ArgumentError slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Iterative(), interpolant_method=Sibson(),
        alpha=0.01, use_cubic_terms=false,
        use_sibson_weight=true,
        tri=tri, z=z)
    itp1 = interpolate(tri, z; derivatives=true)
    itp2 = interpolate(tri, z; derivatives=true)
    ∂11 = differentiate(itp1, 1)
    ∂12 = differentiate(itp1, 2)
    ∂21 = differentiate(itp2, 1)
    ∂22 = differentiate(itp2, 2)
    @test slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22;
        x=5.872, y=3.45, rng=StableRNG(29991),
        method=Iterative(), interpolant_method=Sibson(),
        alpha=0.1, use_cubic_terms=false,
        use_sibson_weight=true,
        tri=tri, z=z)
end

@testset "Check multithreading is working" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(tri, z; generate_derivatives=true)
    ∂1 = differentiate(itp, 1)
    ∂2 = differentiate(itp, 2)
    x = 10rand(100)
    y = 10rand(100)
    @test collect.(∂1(x, y; parallel=true)) ≈ collect.(∂1(x, y; parallel=false)) # not == because of internal rng
    @test collect.(∂1(x, y; interpolant_method = Sibson(1))) ≈ collect.(∂1(x, y; interpolant_method = Sibson(1)))
end