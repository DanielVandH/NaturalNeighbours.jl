using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using ReferenceTests
using CairoMakie
function plot_2d(fig, i, j, title, vals, xg, yg, x, y, show_scatter=true)
    ax = Axis(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    contourf!(ax, xg, yg, reshape(vals, (length(xg), length(yg))), colormap=:viridis, levels=-1:0.05:0, extendlow=:auto, extendhigh=:auto)
    show_scatter && scatter!(ax, vec([(i - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)]), vec([(j - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)]), color=:red, markersize=14)
end
function plot_3d(fig, i, j, title, vals, xg, yg)
    ax = Axis3(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    surface!(ax, xg, yg, reshape(vals, (length(xg), length(yg))), color=vals, colormap=:viridis)
end

@testset "A domain with no holes" begin
    a, b, c, d = 0.0, 1.0, 0.0, 1.0
    nx, ny = 10, 10
    tri = triangulate_rectangle(a, b, c, d, nx, ny)
    f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
    z = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(tri)]
    itp = interpolate(tri, z; derivatives=true)
    xg = LinRange(0, 1, 100)
    yg = LinRange(0, 1, 100)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    exact = f.(_x, _y)
    sibson_vals = itp(_x, _y; method=Sibson())
    triangle_vals = itp(_x, _y; method=Triangle(; allow_cache = false))
    laplace_vals = itp(_x, _y; method=Laplace())
    sibson_1_vals = itp(_x, _y; method=Sibson(1))
    nearest_vals = itp(_x, _y; method=Nearest())
    farin_vals = itp(_x, _y; method=Farin())
    hiyoshi_vals = itp(_x, _y; method=Hiyoshi(2))
    all_vals = (sibson_vals, triangle_vals, laplace_vals, sibson_1_vals, nearest_vals, farin_vals, hiyoshi_vals, exact)
    titles = ("(a): Sibson", "(b): Triangle", "(c): Laplace", "(d): Sibson-1", "(e): Nearest", "(f): Farin", "(g): Hiyoshi", "(h): Exact")
    fig = Figure(fontsize=55)
    for (i, (vals, title)) in enumerate(zip(all_vals, titles))
        plot_2d(fig, 1, i, title, vals, xg, yg, first.(DelaunayTriangulation.each_point(tri)), last.(DelaunayTriangulation.each_point(tri)), !(vals == exact))
        plot_3d(fig, 2, i, " ", vals, xg, yg)
    end
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "../..", "example.png") fig by = psnr_equality(20)
end

@testset "A domain with holes" begin
    R₁ = 2.0
    R₂ = 3.0
    θ = (collect ∘ LinRange)(0, 2π, 250)
    θ[end] = θ[begin]
    x = [
        [R₂ .* cos.(θ)],
        [reverse(R₁ .* cos.(θ))] # inner boundaries are clockwise
    ]
    y = [
        [R₂ .* sin.(θ)],
        [reverse(R₁ .* sin.(θ))] # inner boundaries are clockwise
    ]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes)
    A = get_area(tri)
    D = 6.25e-4
    Tf = (x, y) -> let r = sqrt(x^2 + y^2)
        (R₂^2 - r^2) / (4D) + R₁^2 * log(r / R₂) / (2D)
    end
    _safe_Tf = (x, y) -> let r = sqrt(x^2 + y^2)
        !(R₁ ≤ r ≤ R₂) && return Inf
        return Tf(x, y)
    end
    z = [Tf(x, y) for (x, y) in DelaunayTriangulation.each_point(tri)]
    x = first.(DelaunayTriangulation.each_point(tri))
    y = last.(DelaunayTriangulation.each_point(tri))
    triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
    itp = interpolate(tri, z; derivatives=true)
    xg = LinRange(-R₂, R₂, 75)
    yg = LinRange(-R₂, R₂, 75)
    _x = vec([x for x in xg, _ in yg])
    _y = vec([y for _ in xg, y in yg])
    exact = _safe_Tf.(_x, _y)
    sibson_vals = itp(_x, _y; method=Sibson(), project=false)
    triangle_vals = itp(_x, _y; method=Triangle(), project=false)
    laplace_vals = itp(_x, _y; method=Laplace(), project=false)
    sibson_1_vals = itp(_x, _y; method=Sibson(1), project=false)
    nearest_vals = itp(_x, _y; method=Nearest(), project=false)
    farin_vals = itp(_x, _y; method=Farin(), project=false)
    hiyoshi_vals = itp(_x, _y; method=Hiyoshi(2), project=false)
    all_vals = (sibson_vals, triangle_vals, laplace_vals, sibson_1_vals, nearest_vals, farin_vals, hiyoshi_vals, exact)
    titles = ("(a): Sibson", "(b): Triangle", "(c): Laplace", "(d): Sibson-1", "(e): Nearest", "(f): Farin", "(g): Hiyoshi", "(h): Exact")
    _tri = triangulate([_x'; _y'])
    _triangles = [T[j] for T in each_solid_triangle(_tri), j in 1:3]
    fig = Figure(fontsize=55)
    for (i, (vals, title)) in enumerate(zip(all_vals, titles))
        ax = Axis(fig[1, i], width=600, height=600, title=title, titlealign=:left)
        _vals = copy(vals)
        _vals[isinf.(vals)] .= Inf
        contourf!(ax, _x, _y, _vals, levels=0:50:900)
    end
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "example_constrained.png") fig by=psnr_equality(15)
end