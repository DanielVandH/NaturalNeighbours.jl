using ..NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using ReferenceTests
const NNI = NaturalNeighbours
const DT = DelaunayTriangulation

points = [
    (0.0, 0.0), (-1.0, 1.0), (-0.5, 1.0), (0.0, 1.0), (0.5, 1.0), (1.0, 1.0),
    (1.0, 0.8), (1.0, 0.0), (1.0, -0.5), (1.0, -1.0),
    (0.1, -1.0), (-0.8, -1.0), (-1.0, -1.0),
    (-1.0, -0.7), (-1.0, -0.1), (-1.0, 0.6),
    (-0.1, -0.8), (0.2, -0.8),
    (-0.6, -0.4), (0.9, 0.0), (-0.5, 0.5), (-0.4, 0.6), (-0.1, 0.8)
]
z = zeros(length(points))
z[1] = 1.0
itp = interpolate(points, z, derivatives=true)
vorn = voronoi(itp.triangulation)
xg = LinRange(-1, 1, 250)
yg = LinRange(-1, 1, 250)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
sibson_vals = itp(x, y; method=:sibson) |> x -> reshape(x, (length(xg), length(yg)))
triangle_vals = itp(x, y; method=:triangle) |> x -> reshape(x, (length(xg), length(yg)))
laplace_vals = itp(x, y; method=:laplace) |> x -> reshape(x, (length(xg), length(yg)))
sibson_1_vals = itp(x, y; method=Sibson(1)) |> x -> reshape(x, (length(xg), length(yg)))
nearest = itp(x, y; method=:nearest) |> x -> reshape(x, (length(xg), length(yg)))

function plot_influence(i, j, title, vals, xg, yg, vorn, points)
    ax = Axis(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    contourf!(ax, xg, yg, vals, color=vals, colormap=:viridis, levels=0:0.05:1, extendlow=:auto, extendhigh=:auto)
    voronoiplot!(ax, vorn, strokecolor=:red)
    scatter!(ax, points, color=:red)
    xlims!(ax, -1, 1)
    ylims!(ax, -1, 1)
end
function plot_3d_influence(i, j, title, vals, xg, yg, vorn, points, z)
    ax = Axis3(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    scatter!(ax, first.(points), last.(points), z, color=:red)
    surface!(ax, xg, yg, vals, color=vals, colormap=:viridis, levels=0:0.05:1, extendlow=:auto, extendhigh=:auto)
    xlims!(ax, -1, 1)
    ylims!(ax, -1, 1)
end

fig = Figure(fontsize=36)
plot_influence(1, 1, "(a): Sibson", sibson_vals, xg, yg, vorn, points)
plot_influence(1, 2, "(b): Triangle", triangle_vals, xg, yg, vorn, points)
plot_influence(1, 3, "(c): Laplace", laplace_vals, xg, yg, vorn, points)
plot_influence(1, 4, "(d): Sibson-1", sibson_1_vals, xg, yg, vorn, points)
plot_influence(1, 5, "(e): Nearest", nearest, xg, yg, vorn, points)
plot_3d_influence(2, 1, " ", sibson_vals, xg, yg, vorn, points, z)
plot_3d_influence(2, 2, " ", triangle_vals, xg, yg, vorn, points, z)
plot_3d_influence(2, 3, " ", laplace_vals, xg, yg, vorn, points, z)
plot_3d_influence(2, 4, " ", sibson_1_vals, xg, yg, vorn, points, z)
plot_3d_influence(2, 5, " ", nearest, xg, yg, vorn, points, z)
resize_to_layout!(fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "influence.png") fig
