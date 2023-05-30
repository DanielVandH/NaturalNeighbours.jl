using ..NaturalNeighbours
using DelaunayTriangulation
using StableRNGs
const NNI = NaturalNeighbours
using Random
const DT = DelaunayTriangulation

include(normpath(@__DIR__, "../.", "helper_functions", "point_generator.jl"))

@testset "taylor_neighbourhood" begin
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

@testset "get_λ" begin
    @test NNI.get_λ([1.0, 2.0, 3.0], 2, true) == 2.0
    @test NNI.get_λ([1.0, 2.0, 3.0, 4.0], 5, true) == 1.0
    @test NNI.get_λ([2.3, 5.0], 1, false) == 1.0
    @test NNI.get_λ(1.0, 3, true) == 1.0
end