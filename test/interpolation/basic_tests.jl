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

    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30)
    z = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(tri)]
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

@testset "Sibson(1) errors without gradients" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(tri)]
    itp = interpolate(tri, z; derivatives=false)
    @test_throws ArgumentError("Gradients must be provided for Sibson-1, Farin, or Hiyoshi-2 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true).") itp(0.5, 0.5; method=Sibson(1))
end

@testset "Hiyoshi(2) errors without gradients and Hessians" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(tri)]
    itp = interpolate(tri, z; derivatives=false)
    @test_throws ArgumentError("Gradients and Hessians must be provided for Hiyoshi-2 interpolation. Consider using e.g. interpolate(tri, z; derivatives = true).") itp(0.5, 0.5; method=Hiyoshi(2))
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
        for method in (Sibson(1), Sibson(), Laplace(), Farin(1), Hiyoshi(2), Triangle(), Triangle(; allow_cache = false), Nearest())
            @inferred itp(rand(), rand(); method=method)
            @inferred itp(rand(), rand(); method=method, project=false)
            @inferred itp(rand(Float32), rand(Float32); method=method)
            @inferred itp(rand(Float32), rand(Float32); method=method, project=false)
            @inferred itp(rand(Float32), rand(Float64); method=method)
            @inferred itp(rand(Float32), rand(Float64); method=method, project=false)
            @inferred itp(rand(Float64), rand(Float32); method=method)
            @inferred itp(rand(Float64), rand(Float32); method=method, project=false)
        end
    end
    for method in (Sibson(1), Sibson(), Laplace(), Farin(1), Hiyoshi(2), Triangle(), Nearest())
        p, q = rand(2)
        @test itp1(p, q; method=method) ≈ itp2(p, q; method=method)
        @test itp1(p, q; method=method, project=false) ≈ itp2(p, q; method=method, project=false)
        @test itp1(Float32(p), Float32(q); method=method) ≈ itp2(Float32(p), Float32(q); method=method)
        @test itp1(Float32(p), Float32(q); method=method, project=false) ≈ itp2(Float32(p), Float32(q); method=method, project=false)
        @test itp1(Float32(p), q; method=method) ≈ itp2(Float32(p), q; method=method)
        @test itp1(Float32(p), q; method=method, project=false) ≈ itp2(Float32(p), q; method=method, project=false)
        @test itp1(p, Float32(q); method=method) ≈ itp2(p, Float32(q); method=method)
        @test itp1(p, Float32(q); method=method, project=false) ≈ itp2(p, Float32(q); method=method, project=false)
    end
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
    ch = get_convex_hull_vertices(tri1)
    bad_idx = identify_exterior_points(_itp_xs, _itp_ys, points, ch; tol=1e-3) # boundary effects _really_ matter...
    deleteat!(err, bad_idx)
    @test norm(err) ≈ 0 atol = 1e-2
end


@testset "Test that the derivatives are all zero for missing vertices" begin
    R₁ = 0.2
    R₂ = 1.0
    θ = collect(LinRange(0, 2π, 100))
    θ[end] = 0.0 # get the endpoints to match
    x = [
        [R₂ .* cos.(θ)], # outer first
        [reverse(R₁ .* cos.(θ))] # then inner - reverse to get clockwise orientation
    ]
    y = [
        [R₂ .* sin.(θ)], # 
        [reverse(R₁ .* sin.(θ))]
    ]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes)
    A = get_area(tri)
    refine!(tri; max_area=1e-4A)
    itp = interpolate(tri, ones(DelaunayTriangulation.num_points(tri)), derivatives=true, parallel=false)
    ind = findall(DelaunayTriangulation.each_point_index(tri)) do i
        !DelaunayTriangulation.has_vertex(tri, i)
    end
    for i in ind
        ∇ = NaturalNeighbours.get_gradient(itp, i)
        @test all(iszero, ∇)
        H = NaturalNeighbours.get_hessian(itp, i)
        @test all(iszero, H)
    end
end

@testset "Testing Triangle()'s cache" begin
    tri = triangulate(rand(2, 50))
    method2 = Triangle()
    method2 = Triangle(; allow_cache=false)
    s = Dict{NTuple{3,Int},NTuple{9,Float64}}()
    for T in each_solid_triangle(tri)
        V = DelaunayTriangulation.sort_triangle(T)
        s[V] = NaturalNeighbours._compute_triangle_shape_coefficients(tri, V...)
        @test sum(s[V]) ≈ 1.0
    end
    itp = interpolate(tri, rand(50))
    itp(rand(50), rand(50); method)
    @test isempty(method.s)
    itp(1 / 2, 1 / 2; method=method2)
    @test isempty(method2.s)
    itp(rand(50), rand(50); method=method2)
    @test !isempty(method2.s)
    @test method2.s == s
    tri2 = triangulate(rand(2, 68))
    itp = interpolate(tri2, rand(68))
    itp(rand(50), rand(50); method=method2)
    s = Dict{NTuple{3,Int},NTuple{9,Float64}}()
    for T in each_solid_triangle(tri2)
        V = DelaunayTriangulation.sort_triangle(T)
        s[V] = NaturalNeighbours._compute_triangle_shape_coefficients(tri2, V...)
        @inferred NaturalNeighbours._compute_triangle_shape_coefficients(tri2, V...)
        @test sum(s[V]) ≈ 1.0
    end
    @test method2.s == s # make sure that method2 knows to regenerate the cache
end