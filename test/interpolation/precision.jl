using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs
using LinearAlgebra

include(normpath(@__DIR__, "../.", "helper_functions", "test_functions.jl"))

@testset "Does Sibson-0 reproduce linear functions p ↦ a + bᵀp?" begin
    a = 0.9881
    b = [1.7, 2.3]
    p = [0.3, 0.7]
    f = (x, y) -> a + b' * [x, y]
    xx = LinRange(-10, 10, 25)
    yy = LinRange(-10, 10, 25)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    z = f.(x, y)
    itp = interpolate(x, y, z; derivatives=true)

    xg = LinRange(-10, 10, 250)
    yg = LinRange(-10, 10, 250)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    vals = itp(_x, _y; method=Sibson(0))
    for i in eachindex(vals)
        ξ, η = _x[i], _y[i]
        if DT.distance_to_polygon((ξ, η), get_points(itp.triangulation), get_convex_hull_indices(itp.triangulation)) > 1e-7
            @test vals[i] ≈ f(_x[i], _y[i]) atol = 1e-12
        end
    end
end

@testset "Does Laplace reproduce linear functions p ↦ a + bᵀp?" begin
    a = 0.5673634
    b = [11.7, 62.3]
    p = [0.6, -0.7]
    f = (x, y) -> a + b' * [x, y]
    xx = LinRange(-10, 10, 25)
    yy = LinRange(-10, 10, 25)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    z = f.(x, y)
    itp = interpolate(x, y, z; derivatives=true)

    xg = LinRange(-10, 10, 250)
    yg = LinRange(-10, 10, 250)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    vals = itp(_x, _y; method=Laplace())
    for i in eachindex(vals)
        ξ, η = _x[i], _y[i]
        if DT.distance_to_polygon((ξ, η), get_points(itp.triangulation), get_convex_hull_indices(itp.triangulation)) > 1e-7
            @test vals[i] ≈ f(_x[i], _y[i]) atol = 1e-12
        end
    end
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

@testset "Does Farin reproduce cubics p ↦ a + bᵀx + xᵀQx, Q= [c d; 0 f]?" begin
    a = 0.29912
    b = [1.7, -2.11]
    Q = [2.0 1.01; 0.0 -2.30]
    f = (x, y) -> a + b' * [x, y] + [x, y]' * Q * [x, y]
    xx = LinRange(-10, 10, 25)
    yy = LinRange(-10, 10, 25)
    x = vec([x for x in xx, _ in yy])
    y = vec([y for _ in xx, y in yy])
    z = f.(x, y)
    itp = interpolate(x, y, z; derivatives=true)

    xg = LinRange(-10, 10, 250)
    yg = LinRange(-10, 10, 250)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    vals = itp(_x, _y; method=Farin(1))
    for i in eachindex(vals)
        ξ, η = _x[i], _y[i]
        if DT.distance_to_polygon((ξ, η), get_points(itp.triangulation), get_convex_hull_indices(itp.triangulation)) > 1e-7
            @test vals[i] ≈ f(_x[i], _y[i]) atol = 1e-12
        end
    end
end