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
    @test NNI.get_quadratic_matrix_no_cubic(derivative_cache) == derivative_cache.quadratic_matrix_no_cubic == ElasticMatrix{Float64}(undef, 5, 0)
end

@testset "dwrap" begin
    @test NNI.dwrap(NNI.Direct()) == NNI.Direct()
    @test NNI.dwrap(:direct) == NNI.Direct()
    @test_throws ArgumentError NNI.dwrap(:dir)
    @test NNI.dwrap(NNI.Iterative()) == NNI.Iterative()
    @test NNI.dwrap(:iterative) == NNI.Iterative()
    @test_throws ArgumentError NNI.dwrap(:iter)
end

@testset "show" begin
    tri = triangulate_rectangle(0, 1, 0, 1, 2, 5, add_ghost_triangles=false)
    tri = Triangulation(tri.points, tri.triangles, tri.convex_hull.indices)
    f = (x, y) -> sin(x) + cos(x - y)
    unlock_convex_hull!(tri)
    x = getx.(tri.points)
    y = gety.(tri.points)
    z = f.(x, y)
    ∇ = z .^ 2
    H = z .^ (1 / 5)
    itp = interpolate(tri, z; hessian=∇, gradient=H)
    ∂ = differentiate(itp, 2)
    @test sprint() do io
        Base.show(io, MIME"text/plain"(), ∂)
    end ==
          "Natural Neighbour Differentiator\n    Order: 2\n    z: [1.0, 1.3817732906760363, 0.9689124217106447, 1.5731598536817173, 0.8775825618903728, 1.7190535466982693, 0.7316888688738209, 1.8103834065185413, 0.5403023058681398, 1.8414709848078965]\n    ∇: [1.0, 1.0668106895787572, 0.9937036950756749, 1.094849869590562, 0.9742212470670031, 1.114443051978001, 0.9394318704459826, 1.1260407621773936, 0.8841528765017798, 1.1298817035265263]\n    H: [1.0, 1.909297426825682, 0.9387912809451863, 2.474831925235882, 0.7701511529340699, 2.9551450964158987, 0.5353686008338515, 3.277488078597678, 0.2919265817264289, 3.391015387889364]"
end