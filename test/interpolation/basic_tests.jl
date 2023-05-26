using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs
using LinearAlgebra
using CairoMakie

include(normpath(@__DIR__, "../.", "helper_functions", "test_functions.jl"))

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
    Base.Threads.nthreads() > 1 && @test NNI.get_neighbour_cache(itp, 2) == itp.neighbour_cache[2]
    @test NNI.get_derivative_cache(itp) == itp.derivative_cache
    @test NNI.get_derivative_cache(itp, 1) == itp.derivative_cache[1]
    Base.Threads.nthreads() > 1 && @test NNI.get_derivative_cache(itp, 2) == itp.derivative_cache[2]
    @test NNI.get_gradient(itp) == itp.gradient
    @test !isnothing(NNI.get_gradient(itp))
    @test NNI.get_gradient(itp, 1) == itp.gradient[1]
    @test NNI.get_gradient(itp, 2) == itp.gradient[2]
    @test NNI.get_hessian(itp) == itp.hessian
    @test !isnothing(NNI.get_hessian(itp))
    @test NNI.get_hessian(itp, 1) == itp.hessian[1]
    @test NNI.get_hessian(itp, 2) == itp.hessian[2]
    _itp = interpolate(tri, z; derivatives=false, parallel=false)
    @test NNI.get_gradient(_itp) === nothing
    @test NNI.get_hessian(_itp) === nothing
    @test itp isa NNI.NaturalNeighboursInterpolant
    DT.lock_convex_hull!(tri)
    @test_throws ArgumentError interpolate(tri, z)
    DT.unlock_convex_hull!(tri)
    @test_throws AssertionError interpolate(tri, z[1:end-1])
    w = rand(length(z))
    y = rand(length(z))
    __itp = interpolate(tri, z; derivatives=false, parallel=false, gradient=w)
    @test NNI.get_gradient(__itp) === w
    __itp = interpolate(tri, z; derivatives=false, parallel=false, gradient=w, hessian=y)
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

@testset "Does Sibson-1 reproduce spherical quadratics p ↦ μ(p-a)'(p-a)?" begin
    μ = 0.05
    a = [0.3, 0.7]
    f = (x, y) -> let p = [x, y]
        μ * (p - a)' * (p - a)
    end
    xx = LinRange(a[1] - 5, a[1] + 5, 25)
    yy = LinRange(a[2] - 5, a[2] + 5, 25)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    z = f.(x, y)
    itp = interpolate(x, y, z; derivatives=true)

    xg = LinRange(a[1] - 5, a[1] + 5, 250)
    yg = LinRange(a[2] - 5, a[2] + 5, 250)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    vals = itp(_x, _y; method=Sibson(1))
    for i in eachindex(vals)
        ξ, η = _x[i], _y[i]
        if DT.distance_to_polygon((ξ, η), get_points(itp.triangulation), get_convex_hull_indices(itp.triangulation)) > 1e-7
            @test vals[i] ≈ f(_x[i], _y[i]) atol = 1e-14
        end
    end
end

@testset "Sibson(1) errors without gradients" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(tri, z; derivatives=false)
    @test_throws ArgumentError("Gradients must be provided for Sibson-1 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true).") itp(0.5, 0.5; method=Sibson(1))
end

@testset "Test Float32" begin
    rng = StableRNG(123)
    xs = randn(rng, 100)
    ys = randn(rng, 100)
    tri1 = triangulate([Float32.(xs)'; Float32.(ys)']; rng)
    tri2 = triangulate([xs'; ys']; rng)
    zs = sin.(xs) .* cos.(ys)
    itp1 = interpolate(tri1, Float32.(zs); derivatives=true)
    itp2 = interpolate(tri2, zs; derivatives=true)
    for itp in (itp1, itp2)
        @inferred itp(0.5, 0.5; method=Sibson(1))
        @inferred itp(0.5, 0.5; method=Sibson())
        @inferred itp(0.5, 0.5; method=Laplace())
        @inferred itp(0.5, 0.5; method=Triangle())
        @inferred itp(0.5, 0.5; method=Nearest())
        @inferred itp(0.5, 0.5; method=Sibson(), project=false)
        @inferred itp(0.5, 0.5; method=Laplace(), project=false)
        @inferred itp(0.5, 0.5; method=Triangle(), project=false)
        @inferred itp(0.5, 0.5; method=Nearest(), project=false)
        @inferred itp(0.5, 0.5; method=Sibson(1), project=false)
        @inferred itp(0.5f0, 0.5f0; method=Sibson(1))
        @inferred itp(0.5f0, 0.5f0; method=Sibson())
        @inferred itp(0.5f0, 0.5f0; method=Laplace())
        @inferred itp(0.5f0, 0.5f0; method=Triangle())
        @inferred itp(0.5f0, 0.5f0; method=Nearest())
        @inferred itp(0.5f0, 0.5f0; method=Sibson(), project=false)
        @inferred itp(0.5f0, 0.5f0; method=Laplace(), project=false)
        @inferred itp(0.5f0, 0.5f0; method=Triangle(), project=false)
        @inferred itp(0.5f0, 0.5f0; method=Nearest(), project=false)
        @inferred itp(0.5f0, 0.5f0; method=Sibson(1), project=false)
        @inferred itp(0.5f0, 0.5; method=Sibson(1))
        @inferred itp(0.5f0, 0.5; method=Sibson())
        @inferred itp(0.5f0, 0.5; method=Laplace())
        @inferred itp(0.5f0, 0.5; method=Triangle())
        @inferred itp(0.5f0, 0.5; method=Nearest())
        @inferred itp(0.5f0, 0.5; method=Sibson(), project=false)
        @inferred itp(0.5f0, 0.5; method=Laplace(), project=false)
        @inferred itp(0.5f0, 0.5; method=Triangle(), project=false)
        @inferred itp(0.5f0, 0.5; method=Nearest(), project=false)
        @inferred itp(0.5f0, 0.5; method=Sibson(1), project=false)
        @inferred itp(0.5, 0.5f0; method=Sibson(1))
        @inferred itp(0.5, 0.5f0; method=Sibson())
        @inferred itp(0.5, 0.5f0; method=Laplace())
        @inferred itp(0.5, 0.5f0; method=Triangle())
        @inferred itp(0.5, 0.5f0; method=Nearest())
        @inferred itp(0.5, 0.5f0; method=Sibson(), project=false)
        @inferred itp(0.5, 0.5f0; method=Laplace(), project=false)
        @inferred itp(0.5, 0.5f0; method=Triangle(), project=false)
        @inferred itp(0.5, 0.5f0; method=Nearest(), project=false)
        @inferred itp(0.5, 0.5f0; method=Sibson(1), project=false)
    end
    @test itp1(0.5f0, 0.5f0; method=Sibson(1)) ≈ itp2(0.5f0, 0.5f0; method=Sibson(1))
    @test itp1(0.5f0, 0.5f0; method=Sibson()) ≈ itp2(0.5f0, 0.5f0; method=Sibson())
    @test itp1(0.5f0, 0.5f0; method=Laplace()) ≈ itp2(0.5f0, 0.5f0; method=Laplace())
    @test itp1(0.5f0, 0.5f0; method=Triangle()) ≈ itp2(0.5f0, 0.5f0; method=Triangle())
    @test itp1(0.5f0, 0.5f0; method=Nearest()) ≈ itp2(0.5f0, 0.5f0; method=Nearest())
    @test itp1(0.5f0, 0.5f0; method=Sibson(), project=false) ≈ itp2(0.5f0, 0.5f0; method=Sibson(), project=false)
    @test itp1(0.5f0, 0.5f0; method=Laplace(), project=false) ≈ itp2(0.5f0, 0.5f0; method=Laplace(), project=false)
    @test itp1(0.5f0, 0.5f0; method=Triangle(), project=false) ≈ itp2(0.5f0, 0.5f0; method=Triangle(), project=false)
    @test itp1(0.5f0, 0.5f0; method=Nearest(), project=false) ≈ itp2(0.5f0, 0.5f0; method=Nearest(), project=false)
    test_interpolant(itp1, xs, ys, zs)
    test_interpolant(itp2, xs, ys, zs)

    xrange = LinRange(-3, 3, 1000) .|> Float32
    yrange = LinRange(-3, 3, 1000) .|> Float32
    itp_xs = [xrange[i] for i in 1:length(xrange), j in 1:length(yrange)]
    itp_ys = [yrange[j] for i in 1:length(xrange), j in 1:length(yrange)]
    _itp_xs = vec(itp_xs)
    _itp_ys = vec(itp_ys)
    vals1 = itp1(_itp_xs, _itp_ys; method=Sibson(1))
    vals2 = itp2(_itp_xs, _itp_ys; method=Sibson(1))
    err = abs.(vals1 .- vals2)
    points = get_points(tri1)
    ch = get_convex_hull_indices(tri1)
    bad_idx = identify_exterior_points(_itp_xs, _itp_ys, points, ch; tol=1e-3) # boundary effects _really_ matter...
    deleteat!(err, bad_idx)
    @test norm(err) ≈ 0 atol = 1e-2
end