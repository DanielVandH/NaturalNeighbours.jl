# NaturalNeighbourInterp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/NaturalNeighbourInterp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/NaturalNeighbourInterp.jl/dev/)
[![Build Status](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a package for performing [natural neighbour interpolation](https://en.wikipedia.org/wiki/Natural_neighbor_interpolation) over planar data sets. This method of (scattered data) interpolation takes in some data $X = \{(x_i,y_i)\}\_{i=1}^m \subset \mathbb R^2$ with corresponding data values $Z = \{z_i\}_{i=1}^m$ and constructs a function $f \colon \mathbb R^2 \to \mathbb R$ such that $f(x_i, y_i) = z_i$, $i=1,\ldots,m$, based on the _Voronoi tessellation_ of $X$. We use [DelaunayTriangulation.jl](https://github.com/DanielVandH/DelaunayTriangulation.jl) to construct the Voronoi tessellations. More detail is given in the docs.

# Examples

## Example I: No extrapolation

Let's give some quick examples. The first problem we consider is interpolating $f(x, y) = \sin(xy) - \cos(x-y)\exp[-(x-y)^2]$ for $(x, y) \in [0, 1]^2$. The first step is to define the interpolant, accomplished via `interpolate`.

```julia
using NaturalNeighbourInterp
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
x = vec([(i - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
y = vec([(j - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
z = f.(x, y)
itp = interpolate(x, y, z)
```

This `itp` is our interpolant, and it is callable. To demonstrate this, let's define the points to evaluate the interpolant at.

```julia
xx = LinRange(0, 1, 50)
yy = LinRange(0, 1, 50)
_x = vec([x for x in xx, _ in yy])
_y = vec([y for _ in xx, y in yy])
```

Now we can evaluate the interpolant at these points. While it is possible to do e.g. `itp.(x, y)`, it is more efficient to do `itp(x, y)`. This latter version will use multithreading to handle computing the interpolant at many points, significantly accelerating the computation. Note that doing something like

```julia
Base.Threads.@threads for i in eachindex(x, y)
    z[i] = itp(x[i], y[i]) # DO NOT DO
end
```

will not work as `itp` by itself is not thread-safe. So, with that out of the way, let's evaluate the interpolant. We will use both the Sibson and triangle methods; the Sibson method computes the interpolant weights based on the stolen Voronoi areas, while the triangle method just defines a piecewise linear interpolant over each triangle in the Delaunay triangulation of the point set.

```julia
sibson_vals = itp(_x, _y; method=:sibson) # multithreaded
triangle_vals = itp(_x, _y; method=:triangle)
exact_vals = [f(x, y) for x in xx, y in yy]
sibson_vals = reshape(sibson_vals, (length(xx), length(yy)))
triangle_vals = reshape(triangle_vals, (length(xx), length(yy)))

# Get the errors 
sibson_errs = abs.(sibson_vals .- exact_vals) ./ abs.(exact_vals)
triangle_errs = abs.(triangle_vals .- exact_vals) ./ abs.(exact_vals)
```

Now we can plot the results.

```julia
using CairoMakie
fig = Figure(fontsize=33)
make_ax = (i, j, title) -> begin
    Axis(fig[i, j], title=title, titlealign=:left,
        width=400, height=400,
        xticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]), yticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]))
end
ax1 = make_ax(1, 1, L"(a):$ $ Sibson")
contourf!(ax1, xx, yy, sibson_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax2 = make_ax(1, 2, L"(b):$ $ Triangle")
contourf!(ax2, xx, yy, triangle_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax3 = make_ax(1, 3, L"(c):$ $ Exact")
contourf!(ax3, xx, yy, exact_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax4 = make_ax(2, 3, L"(f):$ $ Data")
tricontourf!(ax4, x, y, z, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax5 = make_ax(2, 1, L"(d):$ $ Sibson error")
contourf!(ax5, xx, yy, sibson_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
ax6 = make_ax(2, 2, L"(e):$ $ Triangle error")
contourf!(ax6, xx, yy, triangle_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
for ax in (ax1, ax2, ax3, ax3, ax4, ax5, ax6)
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)
    scatter!(ax, x, y, markersize=9, color=:red)
end
resize_to_layout!(fig)
fig
```

![Interpolation examples](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/blob/5fabee4777d18117bafe1a55b08ad93994fc1b5a/test/figures/example_1.png)


## Example II: Some extrapolation 

Extrapolation is a difficult problem. We do not currently have the best methods available for this (e.g. with [dynamic ghost points](https://doi.org/10.1016/j.cad.2008.08.007)). Instead, any points that are outside of the convex hull of the boundary are projected onto an edge of the convex hull boundary, or at least the line through that edge, and two-point interpolation is applied to the projected point. Here is an example showing what we can expect from this.

```julia
using NaturalNeighbourInterp, CairoMakie, StableRNGs, DelaunayTriangulation 

## Define the interpolant 
rng = StableRNG(1235)
f = (x, y) -> cos(x^2 / 5) + exp(-(1 / 9) * ((x - y / 9)^2 + (y - x)^2))
x = rand(rng, 33)
y = rand(rng, 33)
z = f.(x, y)
itp = interpolate(x, y, z; rng)

## Points to evaluate at 
xx = LinRange(-1/2, 3/2, 250)
yy = LinRange(-1/2, 3/2, 250)
_x = vec([x for x in xx, _ in yy])
_y = vec([y for _ in xx, y in yy])

## Evaluate the interpolant
sibson_vals = itp(_x, _y; method=:sibson) # multithreaded
triangle_vals = itp(_x, _y; method=:triangle)
exact_vals = [f(x, y) for x in xx, y in yy]
sibson_vals = reshape(sibson_vals, (length(xx), length(yy)))
triangle_vals = reshape(triangle_vals, (length(xx), length(yy)))
sibson_errs = abs.(sibson_vals .- exact_vals) ./ abs.(exact_vals)
triangle_errs = abs.(triangle_vals .- exact_vals) ./ abs.(exact_vals)

## Plot 
fig = Figure(fontsize=33)
make_ax = (i, j, title) -> begin
    Axis(fig[i, j], title=title, titlealign=:left,
        width=400, height=400,
        xticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]), yticks=([0, 0.5, 1], [L"0", L"0.5", L"1"]))
end
ax1 = make_ax(1, 1, L"(a):$ $ Sibson")
contourf!(ax1, xx, yy, sibson_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax2 = make_ax(1, 2, L"(b):$ $ Triangle")
contourf!(ax2, xx, yy, triangle_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax3 = make_ax(1, 3, L"(c):$ $ Exact")
contourf!(ax3, xx, yy, exact_vals, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax4 = make_ax(2, 3, L"(f):$ $ Data")
tricontourf!(ax4, x, y, z, colormap=:viridis, levels=20, colorrange=(-1, 0))
ax5 = make_ax(2, 1, L"(d):$ $ Sibson error")
contourf!(ax5, xx, yy, sibson_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
ax6 = make_ax(2, 2, L"(e):$ $ Triangle error")
contourf!(ax6, xx, yy, triangle_errs, colormap=:viridis, levels=20, colorrange=(0, 0.1))
tri = itp.triangulation
ch_idx = get_convex_hull_indices(tri)
ch_points = [get_point(tri, i) for i in ch_idx]
for ax in (ax1, ax2, ax3, ax3, ax4, ax5, ax6)
    xlims!(ax, -1/2, 3/2)
    ylims!(ax, -1/2, 3/2)
    scatter!(ax, x, y, markersize=9, color=:red)
    lines!(ax, ch_points, color=:black, linewidth=3)
end
resize_to_layout!(fig)
fig
```

![Extrapolation examples](https://github.com/DanielVandH/NaturalNeighbourInterp.jl/blob/5fabee4777d18117bafe1a55b08ad93994fc1b5a/test/figures/example_2.png)

It's not perfect, and it would be nice to eventually have good extrapolation tools. 

## Mathematical Detail

To explain this in more detail, let us make some definitions:

1. _Voronoi tessellation_: Let $X = \{(x_i, y_i)\}_{i=1}^m \subseteq \mathbb R^2$ be a set of points in the plane. We define the _Voronoi tessellation_ $\mathcal V(X) := \{\mathcal V_i\}_{i=1}^m$ to be the partition of $\mathbb R^2$ into convex polygons $\mathcal V_i$ such that all points $(x, y) \in \mathcal V_i$ are closer to the _generator_ $(x_i, y_i)$ than to any other generator $(x_j, y_j)$, $i \neq j$, meaning
$$
\mathcal V_i = \{(x, y) \in \mathbb R^2 \mid (x - x_i)^2 + (y - y_i)^2 \leq (x - x_j)^2 + (y - y_j)^2, i \neq j\},
$$

2. _Natural neighbours_: Two points $(x_i, y_i)$ and $(x_j, y_j)$ are _natural neighbours_ if the intersection $\mathcal V_i \cap \mathcal V_j$ is not empty. The set of natural neighbours to $(x_i, y_i)$ will be denoted $\mathcal N_i$.

Using natural neighbours, we can define a set of generalised barycentric coordinates $\lambda_i(x_0, y_0)$ (also called _weights_ or _local coordinates_) about the natural neighbours $\mathcal N_0$ of some query point $(x_0, y_0)$, allowing us to write:

1. $(x_0, y_0) = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)x_i$. (This is the _local coordinates_ property.)
2. $1 = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)$. (This is the _partition of unity_ property.)
3. For all $i \in \mathcal N_0$, $\lambda_i(x_0, y_0) \geq 0$, meaning $(x_0, y_0)$ is a _convex combination_ of the $x_i \in \mathcal N_i$ or the coordinates $\lambda_i$ are _convex_.

Using these three properties, our interpolation $f$ is defined by 

$$
f(x_0, y_0) = \sum_{i \in \mathcal N_0} \lambda_i(x_0, y_0)z_i.
$$

A basic definition for $\lambda_i(x_0, y_0)$ is through _Sibson's coordinates_, which defines 

$$
\hat\lambda_i = \operatorname{Area}(\mathcal V_0' \cap \mathcal V_i), \quad \lambda_i = \left(\sum_{j \in \mathcal N_0} \hat\lambda_j\right)^{-1}\hat\lambda_i.
$$

Here, $\mathcal V_0'$ is defined as the _virtual tile_ of $(x_0, y_0)$, and is the polygon that would appear in the Voronoi tessellation $\mathcal V(X \cup \{(x_0, y_0)\})$ with $(x_0, y_0)$ inserted. Thus, $\hat\lambda_i$ is the amount of area _stolen_ from $\mathcal V_i$ from the insertion of $(x_0, y_0)$, and $\lambda_i$ is the normalised version.