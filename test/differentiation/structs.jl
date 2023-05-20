using DelaunayTriangulation 
using ..NaturalNeighbours 
const NNI = NaturalNeighbours
const DT = DelaunayTriangulation
using ElasticArrays 

@testset "DerivativeCache" begin
    tri = triangulate_rectangle(0, 10, 0, 10, 101, 101)
    derivative_cache = NNI.DerivativeCache(tri)
    @test NNI.get_iterated_neighbourhood(derivative_cache) == derivative_cache.iterated_neighbourhood == Set{Int64}()
    @test NNI.get_second_iterated_neighbourhood(derivative_cache) == derivative_cache.second_iterated_neighbourhood == Set{Int64}()
    @test NNI.get_linear_matrix(derivative_cache) == derivative_cache.linear_matrix == ElasticMatrix{Float64}(undef, 2, 0)
    @test NNI.get_quadratic_matrix(derivative_cache) == derivative_cache.quadratic_matrix == ElasticMatrix{Float64}(undef, 9, 0)
    @test NNI.get_rhs_vector(derivative_cache) == derivative_cache.rhs_vector == Float64[]
    @test NNI.get_linear_sol(derivative_cache) == derivative_cache.linear_sol == [0.0, 0.0]
    @test NNI.get_quadratic_sol(derivative_cache) == derivative_cache.quadratic_sol == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test NNI.get_quadratic_matrix_no_cubic(derivative_cache) == derivative_cache.quadratic_matrix_no_cubic == ElasticMatrix{Float64}(undef, 5, 0)
    @test NNI.get_quadratic_sol_no_cubic(derivative_cache) == derivative_cache.quadratic_sol_no_cubic == [0.0, 0.0, 0.0, 0.0, 0.0]
end

@testset "dwrap" begin
    @test NNI.dwrap(NNI.Direct()) == NNI.Direct()
    @test NNI.dwrap(:direct) == NNI.Direct()
    @test_throws ArgumentError NNI.dwrap(:dir)
    @test NNI.dwrap(NNI.Iterative()) == NNI.Iterative()
    @test NNI.dwrap(:iterative) == NNI.Iterative()
    @test_throws ArgumentError NNI.dwrap(:iter)
end