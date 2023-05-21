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

@testset "show" begin
    tri = triangulate_rectangle(0, 1, 0, 1, 2, 5, add_ghost_triangles=false)
    tri = Triangulation(tri.points, tri.triangles, tri.convex_hull.indices)
    unlock_convex_hull!(tri)
    x = getx.(tri.points)
    y = gety.(tri.points)
    z = f.(x, y)
    ∇ = z .^ 2
    H = z .^ (1 / 5)
    itp = interpolate(tri, z; hessian=∇, gradient=H)
    @test sprint() do io
        Base.show(io, MIME"text/plain"(), itp)
    end ==
          "Natural Neighbour Interpolant\n    z: [1.0, 1.3817732906760363, 0.9689124217106447, 1.5731598536817173, 0.8775825618903728, 1.7190535466982693, 0.7316888688738209, 1.8103834065185413, 0.5403023058681398, 1.8414709848078965]\n    ∇: [1.0, 1.0668106895787572, 0.9937036950756749, 1.094849869590562, 0.9742212470670031, 1.114443051978001, 0.9394318704459826, 1.1260407621773936, 0.8841528765017798, 1.1298817035265263]\n    H: [1.0, 1.909297426825682, 0.9387912809451863, 2.474831925235882, 0.7701511529340699, 2.9551450964158987, 0.5353686008338515, 3.277488078597678, 0.2919265817264289, 3.391015387889364]"
end