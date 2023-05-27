using ..NaturalNeighbours
using CairoMakie
using ReferenceTests
using StableRNGs

## The data 
rng = StableRNG(123)
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
x = vec([(i - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
y = vec([(j - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
z = f.(x, y)

## The interpolant and grid 
itp = interpolate(x, y, z; derivatives=true)
xg = LinRange(0, 1, 100)
yg = LinRange(0, 1, 100)
_x = vec([x for x in xg, _ in yg])
_y = vec([y for _ in xg, y in yg])
exact = f.(_x, _y)

## Evaluate some interpolants 
sibson_vals = itp(_x, _y; method=Sibson())
triangle_vals = itp(_x, _y; method=Triangle())
laplace_vals = itp(_x, _y; method=Laplace())
sibson_1_vals = itp(_x, _y; method=Sibson(1))
nearest_vals = itp(_x, _y; method=Nearest())
farin_vals = itp(_x, _y; method=Farin())

## Plot 
function plot_2d(fig, i, j, title, vals, xg, yg, x, y, show_scatter=true)
    ax = Axis(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    contourf!(ax, xg, yg, reshape(vals, (length(xg), length(yg))), color=vals, colormap=:viridis, levels=-1:0.05:0, extendlow=:auto, extendhigh=:auto)
    show_scatter && scatter!(ax, x, y, color=:red, markersize=14)
end
function plot_3d(fig, i, j, title, vals, xg, yg)
    ax = Axis3(fig[i, j], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    surface!(ax, xg, yg, reshape(vals, (length(xg), length(yg))), color=vals, colormap=:viridis, levels=-1:0.05:0, extendlow=:auto, extendhigh=:auto)
end

all_vals = (sibson_vals, triangle_vals, laplace_vals, sibson_1_vals, nearest_vals, farin_vals, exact)
titles = ("(a): Sibson", "(b): Triangle", "(c): Laplace", "(d): Sibson-1", "(e): Nearest", "(f): Farin", "(g): Exact")
fig = Figure(fontsize=55)
for (i, (vals, title)) in enumerate(zip(all_vals, titles))
    plot_2d(fig, 1, i, title, vals, xg, yg, x, y, !(vals === exact))
    plot_3d(fig, 2, i, " ", vals, xg, yg)
end
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "example.png") fig
