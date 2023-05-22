# NaturalNeighbours

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/NaturalNeighbours.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/NaturalNeighbours.jl/dev/)
[![Build Status](https://github.com/DanielVandH/NaturalNeighbours.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/NaturalNeighbours.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DanielVandH/NaturalNeighbours.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/NaturalNeighbours.jl)

This is a package for performing [natural neighbour interpolation](https://en.wikipedia.org/wiki/Natural_neighbor_interpolation) over planar data sets (amongst some others, like piecewise linear interpolation via triangles or nearest neighbour interpolation -- see the docs), using [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl) to construct the Voronoi tessellations that represents the spatial information. Most of the work in this package is based on [this great thesis](https://kluedo.ub.rptu.de/frontdoor/deliver/index/docId/2104/file/diss.bobach.natural.neighbor.20090615.pdf).

This is not the only package for interpolation. If the methods available here do not suit your needs, see [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) and the packages it links to in its README.

Here is a quick example of how to use the package, demonstrating the available methods for interpolation. See the docs for more examples, including examples for derivative generation. In this example, note that even though we evaluate the interpolant at $100^2$ points, the runtime is extremely fast thanks to the interpolant being local rather than global.

```julia
using NaturalNeighbours
using CairoMakie
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
nearest = itp(_x, _y; method=Nearest())

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

all_vals = (sibson_vals, triangle_vals, laplace_vals, sibson_1_vals, nearest, exact)
titles = ("(a): Sibson", "(b): Triangle", "(c): Laplace", "(d): Sibson-1", "(e): Nearest", "(f): Exact")
fig = Figure(fontsize=36)
for (i, (vals, title)) in enumerate(zip(all_vals, titles))
    plot_2d(fig, 1, i, title, vals, xg, yg, x, y, !(vals === exact))
    plot_3d(fig, 2, i, " ", vals, xg, yg)
end
resize_to_layout!(fig)
fig

# could keep going and differentiating, etc...
# ∂ = differentiate(itp, 2) -- see the docs.
```

![Interpolation example](https://github.com/DanielVandH/NaturalNeighbours.jl/blob/09708d95ca8b778a84fdf229dc26b51de333b7cd/example.png)
