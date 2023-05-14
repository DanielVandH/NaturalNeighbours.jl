using NaturalNeighbourInterp
using Test

@testset "NaturalNeighbourInterp.jl" begin
    # Write your tests here.
end

using DelaunayTriangulation, Random, CairoMakie, StableRNGs, LinearAlgebra
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

# https://core.ac.uk/reader/36727660
function test_1(x, y)
    return 0.75exp(-((9x - 2)^2 + (9y - 2)^2) / 4) + 0.75exp(-(9x + 1)^2 / 49 - (9y + 1) / 10) + 0.5exp(-((9x - 7)^2 + (9y - 3)^2) / 4) - 0.2exp(-(9x - 4)^2 - (9y - 7)^2)
end
function test_2(x, y)
    return (1 / 9) * (tanh(9y - 9x) + 1)
end
function test_3(x, y)
    return (1.25 + cos(5.4y)) / (6(1 + (3x - 1)^2))
end
function test_4(x, y)
    return (1 / 3) * exp(-(81 / 16) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
end
function test_5(x, y)
    return (1 / 3) * exp(-(81 / 4) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
end
function test_6(x, y)
    return (1 / 9) * (64 - 81 * ((x - 1 / 2)^2 + (y - 1 / 2)^2))^(1 / 2) - 1 / 2
end
function point_set_1()
    A = [0.022703 -0.031021
        0.021701 0.257692
        0.001903 0.494360
        0.039541 0.699342
        0.031583 0.910765
        0.132419 0.050133
        0.125444 0.259297
        0.076758 0.417112
        0.062649 0.655223
        0.095867 0.914652
        0.264560 0.029294
        0.208899 0.266878
        0.171473 0.480174
        0.190921 0.687880
        0.230463 0.904651
        0.366317 0.039695
        0.383239 0.238955
        0.346632 0.490299
        0.387316 0.644523
        0.379536 0.893803
        0.414977 -0.028462
        0.420001 0.226247
        0.485566 0.389142
        0.479258 0.632425
        0.397776 0.848971
        0.053989 0.158674
        0.017513 0.341401
        0.050968 0.578285
        0.048706 0.747019
        0.041878 0.996289
        0.109027 0.091855
        0.093454 0.338159
        0.145187 0.561556
        0.145273 0.752407
        0.069556 0.963242
        0.239164 0.060230
        0.276733 0.369604
        0.226678 0.594059
        0.186765 0.818558
        0.242622 0.980541
        0.385766 0.068448
        0.317909 0.312413
        0.377659 0.519930
        0.381292 0.820379
        0.280351 0.971172
        0.427768 0.156096
        0.466363 0.317509
        0.409203 0.508495
        0.481228 0.751101
        0.402732 0.997873
        0.584869 -0.027195
        0.606389 0.270927
        0.574131 0.425942
        0.599010 0.673378
        0.609697 0.924241
        0.661693 0.025596
        0.639647 0.200834
        0.700118 0.489070
        0.690895 0.669783
        0.671889 0.936610
        0.773694 0.028537
        0.741042 0.193658
        0.730603 0.471423
        0.821453 0.668505
        0.807664 0.847679
        0.842457 0.038050
        0.836692 0.208309
        0.847812 0.433563
        0.917570 0.630738
        0.927987 0.904231
        1.044982 -0.012090
        0.985788 0.269584
        1.012929 0.439605
        1.001985 0.694152
        1.041468 0.868208
        0.573008 0.127243
        0.501389 0.347773
        0.610695 0.608471
        0.538062 0.723524
        0.502619 1.030876
        0.642784 0.070783
        0.670396 0.325984
        0.633359 0.509632
        0.689564 0.775957
        0.683767 1.006451
        0.763533 0.102140
        0.825898 0.323577
        0.808661 0.609159
        0.729064 0.802281
        0.817095 1.051236
        0.868405 0.090205
        0.941846 0.331849
        0.859958 0.591014
        0.859633 0.814484
        0.851280 0.969603
        0.967063 0.133411
        0.967631 0.379528
        0.965704 0.504442
        1.035930 0.745992
        0.947151 0.980141]
    return A[:, 1], A[:, 2]
end
function point_set_2()
    A = [
        0.00 0.00
        0.00 1.00
        0.00 0.50
        0.50 1.00
        0.10 0.15
        0.15 0.30
        0.30 0.35
        0.10 0.75
        0.05 0.45
        1.00 0.00
        1.00 1.00
        0.50 0.00
        1.00 0.50
        0.20 0.10
        0.25 0.20
        0.60 0.25
        0.90 0.35
        0.80 0.40
        0.70 0.20
        0.95 0.90
        0.60 0.65
        0.65 0.70
        0.35 0.85
        0.60 0.85
        0.90 0.80
        0.85 0.25
        0.80 0.65
        0.75 0.85
        0.70 0.90
        0.70 0.65
        0.75 0.10
        0.75 0.35
        0.55 0.95
    ]
    return A[:, 1], A[:, 2]
end
function point_set_3()
    A = [
        0.1375 0.97500
        0.9125 0.98750
        0.7125 0.76250
        0.2250 0.83750
        0.0500 0.41250
        0.4750 0.63750
        0.0500 -0.05000
        0.4500 1.03750
        1.0875 0.55000
        0.5375 0.80000
        0.0375 0.75000
        0.1875 0.57500
        0.7125 0.55000
        0.8500 0.43750
        0.7000 0.31250
        0.2750 0.42500
        0.4500 0.28750
        0.8125 0.18750
        0.4500 -0.03750
        1.0000 0.26250
        0.5000 0.46250
        0.1875 0.26250
        0.5875 0.12500
        1.0500 -0.06125
        0.1000 0.11250
    ]
    return A[:, 1], A[:, 2]
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

@testset "Natural coordinates" begin
    for method in (:sibson, :triangle)
        pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
        tri = triangulate(pts, randomise=false, delete_ghosts=false)
        n = 2500
        pts = random_points_in_convex_hull(tri, n)
        for p in Iterators.flatten((pts, each_point(tri)))
            natural_coordinates = NNI.compute_natural_coordinates(tri, p; method)
            @test sum(NNI.get_coordinates(natural_coordinates)) ≈ 1
            δ = NNI.get_barycentric_deviation(natural_coordinates)
            @test δ ≈ 0 atol = 1e-9
        end

        for _ in 1:50
            pts = [(randn() + rand(), rand() - 2randn()) for _ in 1:250]
            tri = triangulate(pts, delete_ghosts=false)
            n = 5000
            random_points = random_points_in_convex_hull(tri, n)
            for p in Iterators.flatten((random_points, each_point(tri)))
                natural_coordinates = NNI.compute_natural_coordinates(tri, p; method)
                @test sum(NNI.get_coordinates(natural_coordinates)) ≈ 1
                δ = NNI.get_barycentric_deviation(natural_coordinates)
                @test δ ≈ 0 atol = 1e-5
            end
        end
    end
    tri = triangulate(rand(2, 50))
    @test_throws ArgumentError NNI.compute_natural_coordinates(tri, (0.0, 0.0); method=:del)
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

function test_interpolant(itp, x, y, f)
    for method in (:sibson, :triangle)
        for _ in 1:500
            vals = itp(x, y; parallel=false, method)
            vals2 = similar(vals)
            itp(vals2, x, y; parallel=false, method)
            vals3 = itp(x, y; parallel=true, method)
            vals4 = similar(vals3)
            itp(vals4, x, y; parallel=true, method)
            for i in eachindex(x, y)
                _x = x[i]
                _y = y[i]
                _z = f isa Function ? f(_x, _y) : f[i]
                @test all(val -> isapprox(val, _z, rtol=1e-1), (itp(_x, _y; method), vals[i], vals2[i], vals3[i], vals4[i]))
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
    test_interpolant(itp, x, y, f)
    test_interpolant(itp, x, y, z)

    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30, add_ghost_triangles=true)
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(get_points(tri), z)
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
    @test _z ≈ __z
    @test __z ≈ itp(1.8, 0.7; method=:triangle)
    @test __z ≈ itp(1.8, 0.7; method=:sibson)
end

rng = StableRNG(123)
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 15, 15, add_ghost_triangles=true)
z = [f(x, y) for (x, y) in each_point(tri)]
itp = interpolate(get_points(tri), z)
xx = LinRange(0, 1, 50)
yy = LinRange(0, 1, 50)
zz = [itp(x, y) for x in xx, y in yy]
_zz = [f(x, y) for x in xx, y in yy]

x = vec([x for x in xx, _ in yy])
y = vec([y for _ in xx, y in yy])
itp(x, y)
@benchmark $itp($x, $y)


fig = Figure()
ax = Axis(fig[1, 1], aspect=1)
contourf!(ax, xx, yy, _zz, colormap=:viridis, levels=20)
ax = Axis(fig[1, 2], aspect=1)
contourf!(ax, xx, yy, zz, colormap=:viridis, levels=20)
fig

## Example I: No extrapolation 
# Define the interpolant
rng = StableRNG(123)
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
x = vec([(i - 1) / 9 for i in (1, 3, 4,5,8,9,10), j in (1,2,3,5,6,7,9,10)])
y = vec([(j - 1) / 9 for i in (1, 3, 4,5,8,9,10), j in (1,2,3,5,6,7,9,10)])
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
contourf!(ax1, xx, yy, sibson_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax2 = make_ax(1, 2, L"(b):$ $ Triangle")
contourf!(ax2, xx, yy, triangle_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax3 = make_ax(1, 3, L"(c):$ $ Exact")
contourf!(ax3, xx, yy, exact_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax4 = make_ax(2, 3, L"(f):$ $ Data")
tricontourf!(ax4, x, y, z, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax5 = make_ax(2, 1, L"(d):$ $ Sibson error")
contourf!(ax5, xx, yy, sibson_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
ax6 = make_ax(2, 2, L"(e):$ $ Triangle error")
contourf!(ax6, xx, yy, triangle_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
for ax in (ax1, ax2, ax3, ax3, ax4, ax5, ax6)
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)
    scatter!(ax, x, y, markersize=9, color=:red)
end
resize_to_layout!(fig)
fig
save(normpath(@__DIR__, "..", "test", "figures", "example_1.png"), fig)
