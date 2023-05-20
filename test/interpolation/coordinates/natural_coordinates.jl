using ..NaturalNeighbours
using StableRNGs
using DelaunayTriangulation
using Random
using LinearAlgebra
const DT = DelaunayTriangulation
const NNI = NaturalNeighbours

include(normpath(@__DIR__, "../..", "helper_functions", "point_generator.jl"))

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