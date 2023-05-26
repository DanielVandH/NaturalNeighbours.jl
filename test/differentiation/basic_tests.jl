using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Optimization
using OptimizationNLopt
using Random
using StableRNGs
using LinearAlgebra

include(normpath(@__DIR__, "../.", "helper_functions", "slow_derivative_tests.jl"))
include(normpath(@__DIR__, "../.", "helper_functions", "point_generator.jl"))

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
        initial_gradients = NNI.generate_gradients(tri, z, derivative_caches, neighbour_caches; parallel=parallel_derivatives)
        _initial_gradients = deepcopy(initial_gradients)
        for i in eachindex(initial_gradients)
            G1, G2 = estimate_gradient_direct(tri, i, z; use_sibson_weight=true)
            G3 = collect(initial_gradients[i])
            @test G1 ≈ G2 rtol = 1e-4
            @test G1 ≈ G3 rtol = 1e-4
            @test G2 ≈ G3 rtol = 1e-4
        end
        Gpar = NNI.generate_gradients(tri, z, (derivative_caches), neighbour_caches; parallel=true)
        Gser = NNI.generate_gradients(tri, z, (derivative_caches), neighbour_caches; parallel=false)
        @test Gpar == Gser
        ∇par, ℋpar = NNI.generate_derivatives(tri, z, (derivative_caches), neighbour_caches; method, initial_gradients, parallel=true)
        ∇ser, ℋser = NNI.generate_derivatives(tri, z, (derivative_caches), neighbour_caches; method, initial_gradients, parallel=false)
        @test initial_gradients == _initial_gradients # make sure initial_gradients is not modified
        @test ∇par == ∇ser
        @test ℋpar == ℋser
        _α = (0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5)
        flags = zeros(Int64, 8, length(_α))
        for (j, α) in enumerate(_α)
            ∇, ℋ = NNI.generate_derivatives(tri, z, (derivative_caches), neighbour_caches; method, initial_gradients, parallel=true, alpha=α)
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
    itp = interpolate(tri, z; derivatives=true)
    ∂1 = differentiate(itp, 1)
    ∂2 = differentiate(itp, 2)

    x = 10rand(100)
    y = 10rand(100)
    @test collect.(∂1(x, y; parallel=true)) ≈ collect.(∂1(x, y; parallel=false)) rtol = 1e-13 # not == because of internal rng
    @test collect.(∂1(x, y; interpolant_method=Sibson(1), parallel=true)) ≈ collect.(∂1(x, y; parallel=false, interpolant_method=Sibson(1))) rtol = 1e-13
    @test collect.(getindex.(∂2(x, y; parallel=true), 1)) ≈ collect.(getindex.(∂2(x, y; parallel=false), 1)) rtol = 1e-13
    @test collect.(getindex.(∂2(x, y; parallel=true), 2)) ≈ collect.(getindex.(∂2(x, y; parallel=false), 2)) rtol = 1e-13
    @test collect.(getindex.(∂2(x, y; parallel=true, interpolant_method=Sibson(1)), 1)) ≈ collect.(getindex.(∂2(x, y; parallel=false, interpolant_method=Sibson(1)), 1)) rtol = 1e-13
    @test collect.(getindex.(∂2(x, y; parallel=true, interpolant_method=Sibson(1)), 2)) ≈ collect.(getindex.(∂2(x, y; parallel=false, interpolant_method=Sibson(1)), 2)) rtol = 1e-13
end

@testset "Check that Iterative() errors without gradients" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    tri = triangulate(get_points(tri), randomise=false)
    f = (x, y) -> x^2 + y^2 + x^3 * y
    f′ = (x, y) -> [2x + 3x^2 * y; 2y + x^3]
    f′′ = (x, y) -> [2+6x*y 3x^2; 3x^2 2]
    z = [f(x, y) for (x, y) in each_point(tri)]
    itp = interpolate(tri, z; derivatives=false)
    ∂2 = differentiate(itp, 2)
    @test_throws ArgumentError("initial_gradients must be provided for iterative derivative estimation. Consider using e.g. interpolate(tri, z; derivatives = true).") ∂2(0.5, 0.5; method=Iterative())
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
    ∂1 = differentiate(itp1, 2)
    ∂2 = differentiate(itp2, 2)
    for ∂ in (∂1, ∂2)
        @inferred ∂(0.5, 0.5; method=Iterative())
        @inferred ∂(0.5, 0.5; method=Direct())
        @inferred ∂(0.5, 0.5; method=Iterative(), project=false)
        @inferred ∂(0.5, 0.5; method=Direct(), project=false)
        @inferred ∂(0.5f0, 0.5f0; method=Iterative())
        @inferred ∂(0.5f0, 0.5f0; method=Direct())
        @inferred ∂(0.5f0, 0.5f0; method=Iterative(), project=false)
        @inferred ∂(0.5f0, 0.5f0; method=Direct(), project=false)
        @inferred ∂(0.5f0, 0.5; method=Iterative())
        @inferred ∂(0.5f0, 0.5; method=Direct())
        @inferred ∂(0.5f0, 0.5; method=Iterative(), project=false)
        @inferred ∂(0.5f0, 0.5; method=Direct(), project=false)
        @inferred ∂(0.5, 0.5f0; method=Iterative())
        @inferred ∂(0.5, 0.5f0; method=Direct())
        @inferred ∂(0.5, 0.5f0; method=Iterative(), project=false)
        @inferred ∂(0.5, 0.5f0; method=Direct(), project=false)
    end
    let ∇H1 = ∂1(0.5f0, 0.5f0; method=Iterative()), ∇H2 = ∂2(0.5f0, 0.5f0; method=Iterative())
        ∇1, H1 = ∇H1
        ∇2, H2 = ∇H2
        @test collect(∇1) ≈ collect(∇2)
        @test collect(H1) ≈ collect(H2)
    end
    let ∇H1 = ∂1(0.5f0, 0.5f0; method=Direct()), ∇H2 = ∂2(0.5f0, 0.5f0; method=Direct())
        ∇1, H1 = ∇H1
        ∇2, H2 = ∇H2
        @test collect(∇1) ≈ collect(∇2)
        @test collect(H1) ≈ collect(H2)
    end
    let ∇H1 = ∂1(0.5f0, 0.5f0; method=Iterative(), project=false), ∇H2 = ∂2(0.5f0, 0.5f0; method=Iterative(), project=false)
        ∇1, H1 = ∇H1
        ∇2, H2 = ∇H2
        @test collect(∇1) ≈ collect(∇2)
        @test collect(H1) ≈ collect(H2)
    end
    let ∇H1 = ∂1(0.5f0, 0.5f0; method=Direct(), project=false), ∇H2 = ∂2(0.5f0, 0.5f0; method=Direct(), project=false)
        ∇1, H1 = ∇H1
        ∇2, H2 = ∇H2
        @test collect(∇1) ≈ collect(∇2)
        @test collect(H1) ≈ collect(H2)
    end
    ∂1 = differentiate(itp1, 1)
    ∂2 = differentiate(itp2, 1)
    for ∂ in (∂1, ∂2)
        @inferred ∂(0.5, 0.5; method=Iterative())
        @inferred ∂(0.5, 0.5; method=Direct())
        @inferred ∂(0.5, 0.5; method=Iterative(), project=false)
        @inferred ∂(0.5, 0.5; method=Direct(), project=false)
        @inferred ∂(0.5f0, 0.5f0; method=Iterative())
        @inferred ∂(0.5f0, 0.5f0; method=Direct())
        @inferred ∂(0.5f0, 0.5f0; method=Iterative(), project=false)
        @inferred ∂(0.5f0, 0.5f0; method=Direct(), project=false)
        @inferred ∂(0.5f0, 0.5; method=Iterative())
        @inferred ∂(0.5f0, 0.5; method=Direct())
        @inferred ∂(0.5f0, 0.5; method=Iterative(), project=false)
        @inferred ∂(0.5f0, 0.5; method=Direct(), project=false)
        @inferred ∂(0.5, 0.5f0; method=Iterative())
        @inferred ∂(0.5, 0.5f0; method=Direct())
        @inferred ∂(0.5, 0.5f0; method=Iterative(), project=false)
        @inferred ∂(0.5, 0.5f0; method=Direct(), project=false)
    end
    let ∇1 = ∂1(0.5f0, 0.5f0; method=Iterative()), ∇2 = ∂2(0.5f0, 0.5f0; method=Iterative())
        @test collect(∇1) ≈ collect(∇2)
    end
    let ∇1 = ∂1(0.5f0, 0.5f0; method=Direct()), ∇2 = ∂2(0.5f0, 0.5f0; method=Direct())
        @test collect(∇1) ≈ collect(∇2)
    end
    let ∇1 = ∂1(0.5f0, 0.5f0; method=Iterative(), project=false), ∇2 = ∂2(0.5f0, 0.5f0; method=Iterative(), project=false)
        @test collect(∇1) ≈ collect(∇2)
    end
    let ∇1 = ∂1(0.5f0, 0.5f0; method=Direct(), project=false), ∇2 = ∂2(0.5f0, 0.5f0; method=Direct(), project=false)
        @test collect(∇1) ≈ collect(∇2)
    end

    xrange = LinRange(-3, 3, 1000) .|> Float32
    yrange = LinRange(-3, 3, 1000) .|> Float32
    itp_xs = [xrange[i] for i in 1:length(xrange), j in 1:length(yrange)]
    itp_ys = [yrange[j] for i in 1:length(xrange), j in 1:length(yrange)]
    _itp_xs = vec(itp_xs)
    _itp_ys = vec(itp_ys)
    ∂1 = differentiate(itp1, 2)
    ∂2 = differentiate(itp2, 2)

    vals1 = ∂1(_itp_xs, _itp_ys; method=Iterative())
    vals2 = ∂2(_itp_xs, _itp_ys; method=Iterative())

    ∇err = [norm(g1 .- g2) for (g1, g2) in zip(first.(vals1), first.(vals2))]
    Herr = [norm(h1 .- h2) for (h1, h2) in zip(last.(vals1), last.(vals2))]

    points = get_points(tri1)
    ch = get_convex_hull_indices(tri1)
    bad_idx = identify_exterior_points(_itp_xs, _itp_ys, points, ch; tol=1e-2) # boundary effects _really_ matter...
    deleteat!(∇err, bad_idx)
    deleteat!(Herr, bad_idx)

    @test norm(∇err) ≈ 0 atol = 1e-2
    @test norm(Herr) ≈ 0 atol = 1e-1

    ∂1 = differentiate(itp1, 1)
    ∂2 = differentiate(itp2, 1)

    vals1 = ∂1(_itp_xs, _itp_ys; method=Iterative())
    vals2 = ∂2(_itp_xs, _itp_ys; method=Iterative())

    ∇err = [norm(g1 .- g2) for (g1, g2) in zip(first.(vals1), first.(vals2))]

    points = get_points(tri1)
    ch = get_convex_hull_indices(tri1)
    bad_idx = identify_exterior_points(_itp_xs, _itp_ys, points, ch; tol=1e-2) # boundary effects _really_ matter...
    deleteat!(∇err, bad_idx)

    @test norm(∇err) ≈ 0 atol = 1e-1
end