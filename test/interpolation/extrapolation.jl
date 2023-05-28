using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra

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

@testset "Basic extrapolation" begin
    tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 5, 10, add_ghost_triangles=true)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    pts = get_points(tri)
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(pts, z, derivatives=true)

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
    _z = itp(getx(p), gety(p); method=Sibson(1))
    @test _z ≈ itp(getx(p), gety(p); method=Sibson(1))
    @test _z ≈ itp(1.8, 0.7; method=Sibson(1))
    @test _z ≈ itp(1.8, 0.7; method=Farin(1))
    @test _z ≈ itp(1.8, 0.7; method=Hiyoshi(2))

    @test isinf(itp(getx(p), gety(p); method=:triangle, project=false))
    @test isinf(itp(getx(p), gety(p); method=:sibson, project=false))
    @test isinf(itp(getx(p), gety(p); method=:laplace, project=false))
    @test isinf(itp(getx(p), gety(p); method=Sibson(1), project=false))
    @test isinf(itp(getx(p), gety(p); method=Farin(1), project=false))
    @test isinf(itp(getx(p), gety(p); method=Hiyoshi(2), project=false))
    @test isinf(itp(1.8, 0.7; method=:triangle, project=false))
    @test isinf(itp(1.8, 0.7; method=:sibson, project=false))
    @test isinf(itp(1.8, 0.7; method=:laplace, project=false))
    @test isinf(itp(1.8, 0.7; method=Sibson(1), project=false))
    @test isinf(itp(1.8, 0.7; method=Farin(1), project=false))
    @test isinf(itp(1.8, 0.7; method=Hiyoshi(2), project=false))

    ∂ = differentiate(itp, 1)
    @test all(isinf, ∂(getx(p), gety(p); project=false))
    ∂ = differentiate(itp, 2)
    @test all(isinf, ∂(getx(p), gety(p); project=false)[1])
    @test all(isinf, ∂(getx(p), gety(p); project=false)[2])
end