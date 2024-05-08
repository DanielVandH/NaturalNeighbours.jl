using ..NaturalNeighbours
using StableRNGs
using DelaunayTriangulation
using Random
using LinearAlgebra
const DT = DelaunayTriangulation
const NNI = NaturalNeighbours

include(normpath(@__DIR__, "../..", "helper_functions", "point_generator.jl"))

to_mat(H) = [H[1] H[3]; H[3] H[2]]
@testset "Natural coordinates" begin
    for method in (Sibson(), Triangle(), Nearest(), Laplace())
        method = NNI.iwrap(method)
        pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0), (14.0, 8.0), (4.0, 4.0), (10.0, 6.0), (6.0, 2.0), (12.0, 4.0), (0.0, 4.0)]
        tri = triangulate(pts, randomise=false, delete_ghosts=false)
        n = 2500
        pts = random_points_in_convex_hull(tri, n)
        for p in Iterators.flatten((pts, DelaunayTriangulation.each_point(tri)))
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
            for p in Iterators.flatten((random_points, DT.each_point(tri)))
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

function _circular_equality(A, B, by=isequal; kwargs...) # slightly tweaked version of circular_equality from DelaunayTriangulation.jl
    if DT.is_circular(A)
        _A = @views A[begin:(end-1)]
    else
        _A = A
    end
    if DT.is_circular(B)
        _B = @views B[begin:(end-1)]
    else
        _B = B
    end
    same_idx = findmin(abs.(_A[begin] .- _B))[2]
    _mapped_B = circshift(_B, -same_idx + 1)
    return by(_A, _mapped_B; kwargs...)
end

@testset "Test coefficient values for each method" begin # used GeoGebra
    for _ in 1:10 # make sure rng is getting passed consistently 
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
        vorn = voronoi(tri, clip=false)
        q = (5.0, 4.0)
        tri2 = deepcopy(tri)
        add_point!(tri2, q, rng=rng)
        vorn2 = voronoi(tri2, clip=false)
        V = get_polygon(vorn2, DelaunayTriangulation.num_points(tri2))
        AX2 = get_area(vorn2, DelaunayTriangulation.num_points(tri2))

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
        @test _circular_equality(nc.indices, [23, 12, 5, 24, 7, 11])
        @test _circular_equality(nc.coordinates, [K1I1J1, F1G1J1K1, AF1G1D1B1A1, A1B1C1, B1D1E1C1, H1I1E1D1G1] ./ AX, ≈, rtol = 1e-2)

        # Triangle 
        nc = NNI.compute_natural_coordinates(NNI.Triangle(), tri, q; rng=rng)
        V = jump_and_march(tri, q; rng)
        @test _circular_equality(nc.indices, [5, 11, 12])
        @test _circular_equality(nc.coordinates, [0.52, 0.3, 0.18], ≈, rtol = 1e-2)

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
        @test _circular_equality(nc.indices, [23, 12, 5, 24, 7, 11])
        @test _circular_equality(nc.coordinates, [w, ℓ, e, z, g, k], ≈, rtol = 1e-2)

        # Nearest 
        nc = NNI.compute_natural_coordinates(NNI.Nearest(), tri, q; rng=rng)
        @test nc.indices == [5]
        @test nc.coordinates ≈ [1.0]

        # Sibson(1)
        tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
        tri = triangulate(get_points(tri), randomise=false)
        f = (x, y) -> sin(x - y) + cos(x + y)
        z = [f(x, y) for (x, y) in DT.each_point(tri)]
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

        # Farin(1)
        tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
        tri = triangulate(get_points(tri), randomise=false)
        f = (x, y) -> sin(x - y) + cos(x + y)
        z = [f(x, y) for (x, y) in DT.each_point(tri)]
        itp = interpolate(tri, z; derivatives=true)
        q = (5.37841, 1.3881)
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q)
        ∇ = NNI.get_gradient(itp)
        λ = NNI.get_coordinates(nc)
        N₀ = NNI.get_indices(nc)
        @test NNI.is_bezier_point(1, 1, 1)
        @test !NNI.is_bezier_point(1, 2, 1)
        @test NNI.is_bezier_edge(1, 2, 2)
        @test !NNI.is_bezier_edge(1, 2, 3)
        @test NNI.is_bezier_face(1, 2, 3)
        @test !NNI.is_bezier_face(1, 2, 2)
        @test NNI.find_bezier_edge(1, 1, 2) == (1, 2)
        @test NNI.find_bezier_edge(1, 2, 1) == (1, 2)
        @test NNI.find_bezier_edge(1, 2, 2) == (2, 1)
        @test collect(NNI.bezier_point_contribution(1, N₀, z)) ≈ [z[N₀[1]], 6]
        @test collect(NNI.bezier_edge_contribution(tri, 2, 5, N₀, ∇, z)) ≈ [z[N₀[2]] + (1 / 3) * dot(get_point(tri, N₀[5]) .- get_point(tri, N₀[2]), ∇[N₀[2]]), 2]
        @test NNI.bezier_face_contribution(tri, 2, 3, 6, N₀, ∇, z)[1] ≈ (1 / 3) * (
            z[N₀[2]] +
            z[N₀[3]] +
            z[N₀[6]]
        ) + (1 / 12) * (
            dot(get_point(tri, N₀[3]) .- get_point(tri, N₀[2]), ∇[N₀[2]]) +
            dot(get_point(tri, N₀[6]) .- get_point(tri, N₀[2]), ∇[N₀[2]]) +
            dot(get_point(tri, N₀[2]) .- get_point(tri, N₀[3]), ∇[N₀[3]]) +
            dot(get_point(tri, N₀[6]) .- get_point(tri, N₀[3]), ∇[N₀[3]]) +
            dot(get_point(tri, N₀[2]) .- get_point(tri, N₀[6]), ∇[N₀[6]]) +
            dot(get_point(tri, N₀[3]) .- get_point(tri, N₀[6]), ∇[N₀[6]])
        )
        @test NNI.bezier_face_contribution(tri, 2, 3, 6, N₀, ∇, z)[2] == 1
        @test collect(NNI.get_contrib(tri, 1, 1, 1, N₀, ∇, z)) ≈ collect(NNI.bezier_point_contribution(1, N₀, z))
        @test collect(NNI.get_contrib(tri, 1, 2, 2, N₀, ∇, z)) ≈ collect(NNI.bezier_edge_contribution(tri, 2, 1, N₀, ∇, z))
        @test collect(NNI.get_contrib(tri, 1, 2, 3, N₀, ∇, z)) ≈ collect(NNI.bezier_face_contribution(tri, 1, 2, 3, N₀, ∇, z))
        n = length(λ)
        s1 = 0.0
        for i in 1:n, j in 1:n, k in 1:n
            s1 += NNI.get_contrib(tri, i, j, k, N₀, ∇, z)[1] * λ[i] * λ[j] * λ[k]
        end
        @test s1 ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test NNI._eval_interp(Farin(1), itp, q, NNI.NaturalNeighboursCache(tri)) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test NNI._eval_natural_coordinates(Farin(1), nc, z, ∇, tri) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test itp(q..., method=Farin(1)) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test itp(5.0, 5.0, method=Farin(1)) ≈ f(5.0, 5.0)
        @test itp(5.0632, 5.0632, method=Farin(1)) ≈ f(5.0632, 5.0632) rtol = 1e-3

        tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 300, 300)
        z = [f(x, y) for (x, y) in DT.each_point(tri)]
        xx = LinRange(0, 1, 50)
        yy = LinRange(0, 1, 50)
        x = vec([x for x in xx, _ in yy])
        y = vec([y for _ in xx, y in yy])
        itp = interpolate(get_points(tri), z; derivatives=true)
        q = (x[52], y[52])
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q)
        ∇ = NNI.get_gradient(itp)
        λ = NNI.get_coordinates(nc)
        N₀ = NNI.get_indices(nc)
        n = length(λ)
        s1 = 0.0
        for i in 1:n, j in 1:n, k in 1:n
            s1 += NNI.get_contrib(tri, i, j, k, N₀, ∇, z)[1] * λ[i] * λ[j] * λ[k]
        end
        @test s1 ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇) rtol = 1e-5
        @test NNI._eval_interp(Farin(1), itp, q, NNI.NaturalNeighboursCache(tri)) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test NNI._eval_natural_coordinates(Farin(1), nc, z, ∇, tri) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test itp(q..., method=Farin(1)) ≈ NNI._compute_farin_coordinates(nc, tri, z, ∇)
        @test itp(q..., method=Farin(1)) ≈ f(q...) rtol = 1e-4

        # Hiyoshi(2)
        tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
        tri = triangulate(get_points(tri), randomise=false)
        f = (x, y) -> sin(x - y) + cos(x + y)
        z = [f(x, y) for (x, y) in DT.each_point(tri)]
        itp = interpolate(tri, z; derivatives=true)
        q = (5.37841, 1.3881)
        nc = NNI.compute_natural_coordinates(NNI.Sibson(), tri, q)
        ∇ = NNI.get_gradient(itp)
        λ = NNI.get_coordinates(nc)
        N₀ = NNI.get_indices(nc)
        H = NNI.get_hessian(itp)
        _z(i) = z[N₀[i]]
        _z(i, j) = dot(∇[N₀[i]], get_point(tri, N₀[j]) .- get_point(tri, N₀[i]))
        _z(i, j, k) = collect(get_point(tri, N₀[j]) .- get_point(tri, N₀[i]))' * to_mat(H[N₀[i]]) * collect(get_point(tri, N₀[k]) .- get_point(tri, N₀[i]))
        @test NNI._hiyoshi_case_1(1, N₀, z) ≈ _z(1)
        let i = 1, j = 2
            @test NNI._hiyoshi_case_2(tri, i, j, N₀, ∇, z) ≈ _z(i) + _z(i, j) / 5
        end
        let i = 1, j = 2
            @test NNI._hiyoshi_case_3(tri, i, j, N₀, ∇, H, z) ≈ _z(i) + 2_z(i, j) / 5 + _z(i, j, j) / 20
        end
        let i = 1, j = 2, k = 3
            @test NNI._hiyoshi_case_4(tri, i, j, k, N₀, ∇, H, z) ≈ _z(i) + (_z(i, j) + _z(i, k)) / 5 + _z(i, j, k) / 20
        end
        let i = 1, j = 2, k = 3
            @test NNI._hiyoshi_case_5(tri, i, j, k, N₀, ∇, H, z) ≈
                  (13 / 30) * (_z(i) + _z(j)) + (2 / 15) * _z(k) +
                  (1 / 9) * (_z(i, j) + _z(j, i)) + (7 / 90) * (_z(i, k) + _z(j, k)) +
                  (2 / 45) * (_z(k, i) + _z(k, j)) + (1 / 45) * (_z(i, j, k) + _z(j, i, k) + _z(k, i, j))
        end
        let i = 1, j = 2, k = 3, ℓ = 4
            @test NNI._hiyoshi_case_6(tri, i, j, k, ℓ, N₀, ∇, H, z) ≈
                  (1 / 2) * _z(i) + (1 / 6) * (_z(j) + _z(k) + _z(ℓ)) +
                  (7 / 90) * (_z(i, j) + _z(i, k) + _z(i, ℓ)) +
                  (2 / 45) * (_z(j, i) + _z(k, i) + _z(ℓ, i)) +
                  (1 / 30) * (_z(j, k) + _z(j, ℓ) + _z(k, j) + _z(k, ℓ) + _z(ℓ, j) + _z(ℓ, k)) +
                  (1 / 90) * (_z(i, j, k) + _z(i, j, ℓ) + _z(i, k, ℓ)) +
                  (1 / 90) * (_z(j, i, k) + _z(j, i, ℓ) + _z(k, i, j) + _z(k, i, ℓ) + _z(ℓ, i, j) + _z(ℓ, i, k)) +
                  (1 / 180) * (_z(j, k, ℓ) + _z(k, j, ℓ) + _z(ℓ, j, k))
        end
        let i = 1, j = 2, k = 3, ℓ = 4, m = 5
            @test NNI._hiyoshi_case_7(tri, i, j, k, ℓ, m, N₀, ∇, H, z) ≈
                  (_z(i) + _z(j) + _z(k) + _z(ℓ) + _z(m)) / 5 +
                  (1 / 30) * (_z(i, j) + _z(i, k) + _z(i, ℓ) + _z(i, m) + _z(j, i) + _z(j, k) + _z(j, ℓ) +
                              _z(j, m) + _z(k, i) + _z(k, j) + _z(k, ℓ) + _z(k, m) + _z(ℓ, i) + _z(ℓ, j) + _z(ℓ, k) + _z(ℓ, m) +
                              _z(m, i) + _z(m, j) + _z(m, k) + _z(m, ℓ)) +
                  (1 / 180) * (_z(i, j, k) + _z(i, j, ℓ) + _z(i, j, m) + _z(i, k, ℓ) + _z(i, k, m) + _z(i, ℓ, m) +
                               _z(j, i, ℓ) + _z(j, i, k) + _z(i, i, m) + _z(j, k, ℓ) + _z(j, k, m) + _z(j, ℓ, m) + _z(k, i, j) +
                               _z(k, i, ℓ) + _z(k, i, m) + _z(k, j, ℓ) + _z(k, j, m) + _z(k, ℓ, m) + _z(ℓ, i, j) + _z(ℓ, i, k) +
                               _z(ℓ, i, m) + _z(ℓ, j, k) + _z(ℓ, j, m) + _z(ℓ, k, m) + _z(m, i, j) + _z(m, i, k) + _z(m, i, ℓ) +
                               _z(m, j, k) + _z(m, j, ℓ) + _z(m, k, ℓ))
        end
        n = length(λ)
        global ss = 0.0
        for a in 1:n, b in 1:n, c in 1:n, d in 1:n, e in 1:n
            prod = λ[a] * λ[b] * λ[c] * λ[d] * λ[e]
            (i, j, k, ℓ, m), case = NNI.group_sort(a, b, c, d, e)
            local s = 0.0
            if case == 1
                let i = i
                    s += _z(i)
                    s *= prod
                end
            elseif case == 2
                let i = ℓ, j = m
                    s += _z(i) + _z(i, j) / 5
                    s *= prod
                end
            elseif case == 3
                let i = i, j = m
                    s += _z(i) + 2_z(i, j) / 5 + _z(i, j, j) / 20
                    s *= prod
                end
            elseif case == 4
                let i = i, j = ℓ, k = m
                    s += _z(i) + (_z(i, j) + _z(i, k)) / 5 + _z(i, j, k) / 20
                    s *= prod
                end
            elseif case == 5
                let i = i, j = k, k = m
                    s += (13 / 30) * (_z(i) + _z(j)) + (2 / 15) * _z(k) +
                         (1 / 9) * (_z(i, j) + _z(j, i)) + (7 / 90) * (_z(i, k) + _z(j, k)) +
                         (2 / 45) * (_z(k, i) + _z(k, j)) + (1 / 45) * (_z(i, j, k) + _z(j, i, k) + _z(k, i, j))
                    s *= prod
                end
            elseif case == 6
                let i = i, j = k, k = ℓ, ℓ = m
                    s += (1 / 2) * _z(i) + (1 / 6) * (_z(j) + _z(k) + _z(ℓ)) +
                         (7 / 90) * (_z(i, j) + _z(i, k) + _z(i, ℓ)) +
                         (2 / 45) * (_z(j, i) + _z(k, i) + _z(ℓ, i)) +
                         (1 / 30) * (_z(j, k) + _z(j, ℓ) + _z(k, j) + _z(k, ℓ) + _z(ℓ, j) + _z(ℓ, k)) +
                         (1 / 90) * (_z(i, j, k) + _z(i, j, ℓ) + _z(i, k, ℓ)) +
                         (1 / 90) * (_z(j, i, k) + _z(j, i, ℓ) + _z(k, i, j) + _z(k, i, ℓ) + _z(ℓ, i, j) + _z(ℓ, i, k)) +
                         (1 / 180) * (_z(j, k, ℓ) + _z(k, j, ℓ) + _z(ℓ, j, k))
                    s *= prod
                end
            elseif case == 7
                s += (_z(i) + _z(j) + _z(k) + _z(ℓ) + _z(m)) / 5 +
                     (1 / 30) * (_z(i, j) + _z(i, k) + _z(i, ℓ) + _z(i, m) + _z(j, i) + _z(j, k) + _z(j, ℓ) +
                                 _z(j, m) + _z(k, i) + _z(k, j) + _z(k, ℓ) + _z(k, m) + _z(ℓ, i) + _z(ℓ, j) + _z(ℓ, k) + _z(ℓ, m) +
                                 _z(m, i) + _z(m, j) + _z(m, k) + _z(m, ℓ)) +
                     (1 / 180) * (_z(i, j, k) + _z(i, j, ℓ) + _z(i, j, m) + _z(i, k, ℓ) + _z(i, k, m) + _z(i, ℓ, m) +
                                  _z(j, i, ℓ) + _z(j, i, k) + _z(i, i, m) + _z(j, k, ℓ) + _z(j, k, m) + _z(j, ℓ, m) + _z(k, i, j) +
                                  _z(k, i, ℓ) + _z(k, i, m) + _z(k, j, ℓ) + _z(k, j, m) + _z(k, ℓ, m) + _z(ℓ, i, j) + _z(ℓ, i, k) +
                                  _z(ℓ, i, m) + _z(ℓ, j, k) + _z(ℓ, j, m) + _z(ℓ, k, m) + _z(m, i, j) + _z(m, i, k) + _z(m, i, ℓ) +
                                  _z(m, j, k) + _z(m, j, ℓ) + _z(m, k, ℓ))
                s *= prod
            end
            ss += s
        end
        @test ss ≈ NNI._compute_hiyoshi_coordinates(nc, tri, z, ∇, H) rtol = 1e-3
        @test NNI._eval_interp(Hiyoshi(2), itp, q, NNI.NaturalNeighboursCache(tri)) ≈ NNI._compute_hiyoshi_coordinates(nc, tri, z, ∇, H)
        @test NNI._eval_natural_coordinates(Hiyoshi(2), nc, z, ∇, H, tri) ≈ NNI._compute_hiyoshi_coordinates(nc, tri, z, ∇, H)
        @test itp(q..., method=Hiyoshi(2)) ≈ NNI._compute_hiyoshi_coordinates(nc, tri, z, ∇, H)
        @test itp(q..., method=Hiyoshi(2)) ≈ f(q...) rtol = 1e-4
    end
end

@testset "show" begin
    interpolation_point = (0.3, 0.7)
    indices = [1, 4, 6, 10]
    coordinates = [0.371, 0.392, 0.4991, 491.20]
    triangulation = triangulate(rand(2, 50))
    nc = NNI.NaturalCoordinates(coordinates, indices, interpolation_point, triangulation)
    @test sprint() do io
        Base.show(io, MIME"text/plain"(), nc)
    end == "NaturalCoordinates{Float64,Int64}\n    u: (0.3, 0.7)\n    λ: [0.371, 0.392, 0.4991, 491.2]\n    k: [1, 4, 6, 10]"
end