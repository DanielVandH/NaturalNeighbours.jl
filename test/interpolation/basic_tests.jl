using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs

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