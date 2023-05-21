using ..NaturalNeighbours
using DelaunayTriangulation
const NNI = NaturalNeighbours

@testset "iwrap" begin
    @test NNI.iwrap(NNI.Sibson()) == NNI.Sibson()
    @test NNI.iwrap(NNI.Triangle()) == NNI.Triangle()
    @test NNI.iwrap(NNI.Nearest()) == NNI.Nearest()
    @test NNI.iwrap(NNI.Laplace()) == NNI.Laplace()
    @test NNI.iwrap(:sibson) == NNI.Sibson()
    @test NNI.iwrap(:triangle) == NNI.Triangle()
    @test NNI.iwrap(:nearest) == NNI.Nearest()
    @test NNI.iwrap(:laplace) == NNI.Laplace()
    @test_throws ArgumentError NNI.iwrap(:lap)

    @test NNI.iwrap(NNI.Sibson(1)) == NNI.Sibson(1)
    @test NNI.iwrap(NNI.Triangle(1)) == NNI.Triangle(0)
    @test_throws ArgumentError NNI.Sibson(5)
    @test NNI.iwrap(NNI.Laplace(1)) == NNI.Laplace(0)
end