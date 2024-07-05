```@meta
CurrentModule = NaturalNeighbours
```

# Differentiation Example 

The purpose of this example is to explore derivative generation. For this, it is important to note that we are thinking of _generating_ derivatives rather than _estimating_ them: Following [Alfeld (1989)](https://doi.org/10.1016/B978-0-12-460515-2.50005-6), derivative generation only seeks to find derivatives that best fit our assumptions of the data, i.e. that give a most satisfactory interpolant, rather than trying to find exact derivative values. The complete quote for this by [Alfeld (1989)](https://doi.org/10.1016/B978-0-12-460515-2.50005-6) is below:

> It seems inevitable that in order to obtain an interpolant that is both local and smooth one has to supply derivative data. Typically, such data are not part of the interpolation problem and have to be made up from existing functional data. This process is usually referred as derivative estimation, but this is probably a misnomer. The objective is not to estimate existing but unknown values of derivatives. Instead, it is to generate values that will yield a satisfactory interpolant. Even if an underlying primitive function did exist it might be preferable to use derivative values that differ from the exact ones. (For example, a maximum error might be decreased by using the "wrong" derivative values.) Therefore, I prefer the term derivative generation rather than derivative estimation.

For the purpose of this exploration, we use Franke's test function. This function, introduced by [Franke and Nielson (1980)](https://doi.org/10.1002/nme.1620151110), is given by 

```math
\begin{align*}
f(x, y) &= \frac34\exp\left\{-\frac{(9x-2)^2 + (9y-2)^2}{4}\right\} + \frac34\exp\left\{-\frac{(9x+1)^2}{49} - \frac{9y+1}{10}\right\} \\
&+ \frac12\exp\left\{-\frac{(9x-7)^2 + (9y-3)^2}{4}\right\} - \frac15\exp\left\{-(9x-4)^2 - (9y-7)^2\right\}.
\end{align*}
```

This function, and its derivatives, are defined below.

```julia
f = (x, y) -> 0.75 * exp(-((9 * x - 2)^2 + (9 * y - 2)^2) / 4) + 0.75 * exp(-(9 * x + 1)^2 / 49 - (9 * y + 1) / 10) + 0.5 * exp(-((9 * x - 7)^2 + (9 * y - 3)^2) / 4) - 0.2 * exp(-(9 * x - 4)^2 - (9 * y - 7)^2)
f′ = (x, y) -> [(exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * x - 72)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * x) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * x) / 2 - 63 / 2)) / 2 - (3 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10) * ((162 * x) / 49 + 18 / 49)) / 4
    (exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * y - 126)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * y) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * y) / 2 - 27 / 2)) / 2 - (27 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)) / 40]
f′′ = (x, y) -> [(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/98-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49)^2)/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)^2)/5 (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5
    (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5 (243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/400+(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*y)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*y)/2-27/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*y-126)^2)/5]
```

Here is the surface of $f$ along with its derivatives.

```julia
using CairoMakie
function plot_f(fig, x, y, vals, title, i, show_3d=true, zlabel="z")
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = contourf!(ax, x, y, vals, colormap=:viridis,  extendhigh=:auto)
    if show_3d
        ax = Axis3(fig[2, i], xlabel="x", ylabel="y", zlabel=zlabel, width=600, height=600, title=" ", titlealign=:left, azimuth=0.49)
        surface!(ax, x, y, vals, color=vals, colormap=:viridis)
    end
    return c
end

x = LinRange(0, 1, 100)
y = LinRange(0, 1, 100)
z = [f(x, y) for x in x, y in y]
∇ = [f′(x, y) for x in x, y in y]
∇₁ = first.(∇)
∇₂ = last.(∇)
H = [f′′(x, y) for x in x, y in y]
H₁₁ = getindex.(H, 1)
H₁₂ = getindex.(H, 2)
H₂₂ = getindex.(H, 4)

fig = Figure(fontsize = 36)
plot_f(fig, x, y, z, "(a): f", 1, true, "z")
plot_f(fig, x, y, ∇₁, "(b): ∂f/∂x", 2, true, "∂f/∂x")
plot_f(fig, x, y, ∇₂, "(c): ∂f/∂y", 3, true,  "∂f/∂y")
plot_f(fig, x, y, H₁₁, "(d): ∂²f/∂x²", 4, true, "∂²f/∂x²")
plot_f(fig, x, y, H₂₂, "(f): ∂²f/∂y²", 5, true, "∂²f/∂y²")
plot_f(fig, x, y, H₁₂, "(e): ∂²f/∂x∂y", 6, true,  "∂²f/∂x∂y")
resize_to_layout!(fig)
fig
```

```@raw html
<figure>
    <img src='../figures/differentiation_exact_surfaces.png', alt'Plots of the interpolants'><br>
</figure>
```

For our analysis, we use the following data set:

```julia
using StableRNGs 
using DelaunayTriangulation 
using CairoMakie
rng = StableRNG(9199)
x = rand(rng, 500)
y = rand(rng, 500)
z = f.(x, y)
tri = triangulate([x'; y'])
vorn = voronoi(tri)

fig = Figure(fontsize=36, size=(1800, 600))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", width=600, height=600, title="(a): Data and triangulation", titlealign=:left)
scatter!(ax, x, y, color=:black, markersize=9)
triplot!(ax, tri, strokecolor=:black, strokewidth=2, show_convex_hull=false)
voronoiplot!(ax, vorn, strokecolor=:blue)
xlims!(ax, 0, 1)
ylims!(ax, 0, 1)

ax = Axis3(fig[1, 2], xlabel="x", ylabel="y", zlabel="z", width=600, height=600, azimuth=0.25, title="(b): Function values", titlealign=:left)
triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
surface!(ax, x, y, z)
scatter!(ax, x, y, z, color=:black, markersize=9)

colgap!(fig.layout, 1, 75)
resize_to_layout!(fig)
fig
```

```@raw html
<figure>
    <img src='../figures/example_data.png', alt'Plots of the data'><br>
</figure>
```

# Generation at the Data Sites 

To start with the example, we consider generating the derivatives at the data sites.

## Gradients 

Let us first estimate some gradients using the direct method.

```julia
using NaturalNeighbours
using DelaunayTriangulation 
using CairoMakie 
using LinearAlgebra

function plot_f2(fig, x, y, vals, title, i, tri, levels, show_3d=true, zlabel="z")
    triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = tricontourf!(ax, x, y, vals, triangulation=triangles', colormap=:viridis, extendhigh=:auto, levels=levels)
    if show_3d
        ax = Axis3(fig[2, i], xlabel="x", ylabel="y", zlabel=zlabel, width=600, height=600, title=" ", titlealign=:left, azimuth=0.49)
        mesh!(ax, hcat(x, y, vals), triangles, color=vals, colormap=:viridis, colorrange=extrema(levels))
    end
    return c
end
function plot_gradients(∇g, tri, f′, x, y)
    ∇g1 = first.(∇g)
    ∇g2 = last.(∇g)
    fig = Figure(fontsize=36, size=(2400, 600))
    plot_f2(fig, x, y, ∇g1, "(a): ∂f̂/∂x", 1, tri, -3.5:0.5:3.0, true, "∂f̂/∂x")
    plot_f2(fig, x, y, ∇g2, "(b): ∂f̂/∂y", 3, tri, -3.5:0.5:3.0, true, "∂f̂/∂y")
    plot_f2(fig, x, y, getindex.(f′.(x, y), 1), "(c): ∂f/∂x", 2, tri, -3.5:0.5:3.0, true, "∂f/∂x")
    plot_f2(fig, x, y, getindex.(f′.(x, y), 2), "(d): ∂f/∂y", 4, tri, -3.5:0.5:3.0, true, "∂f/∂y")
    plot_f2(fig, x, y, norm.(collect.(∇g) .- f′.(x, y)), "(e): Gradient error", 5, tri, 0:0.1:0.5, true, "|∇ε|")
    resize_to_layout!(fig)
    ε = 100sqrt(sum((norm.(collect.(∇g) .- f′.(x, y))) .^ 2) / sum(norm.(∇g) .^ 2))
    return fig, ε
end

points = [x'; y']
z = f.(x, y)
tri = triangulate(points)
∇g = generate_gradients(tri, z)
fig, ε = plot_gradients(∇g, tri, f′, x, y)
```

```julia-repl
julia> ε
10.251180094083372
```

```@raw html
<figure>
    <img src='../figures/gradient_data.png', alt'Gradients'><br>
</figure>
```

A 10% error is not terrible, and the derivatives we obtain are reasonable.

Let's also look at the results when we jointly estimate the gradients and Hessians (this is the default option).

```julia
∇gr, _ = generate_derivatives(tri, z; method=Direct())
fig, ε = plot_gradients(∇gr, tri, f′, x, y)
```

```julia-repl
julia> ε
7.717791597731752
```

```@raw html
<figure>
    <img src='../figures/joint_gradient_data.png', alt'Joint Gradients'><br>
</figure>
```

The figures are smoother, and the error has now decreased to 7.7% -- an improvement. We could also try and see what happens if we use the `Iterative()` approach, using the first gradients we got as the initial gradients, or we could see what happens with `use_cubic_terms=false` for the `Direct()` method, but we won't show that here.

## Hessians 

Let's now look at estimating Hessians. We first consider the direct approach, including cubic terms.

```julia
to_mat(H::NTuple{3,Float64}) = [H[1] H[3]; H[3] H[2]]
function plot_hessians(H, tri, f′′, x, y)
    H₁₁ = getindex.(H, 1)
    H₁₂ = getindex.(H, 3)
    H₂₂ = getindex.(H, 2)

    fig = Figure(fontsize=36, size=(2400, 600))
    plot_f2(fig, x, y, H₁₁, "(a): ∂²f̂/∂x²", 1, tri, -35:5:30, true, "∂²f̂/∂x²")
    plot_f2(fig, x, y, H₂₂, "(c): ∂²f̂/∂y²", 3, tri, -35:5:30, true, "∂²f̂/∂y²")
    plot_f2(fig, x, y, H₁₂, "(e): ∂²f̂/∂x∂y", 5, tri, -35:5:30, true, "∂²f̂/∂x∂y")
    plot_f2(fig, x, y, getindex.(f′′.(x, y), 1), "(b): ∂²f/∂x²", 2, tri, -35:5:30, true, "∂²f/∂x²")
    plot_f2(fig, x, y, getindex.(f′′.(x, y), 4), "(d): ∂²f/∂y²", 4, tri, -35:5:30, true, "∂²f/∂y²")
    plot_f2(fig, x, y, getindex.(f′′.(x, y), 2), "(f): ∂²f/∂x∂y", 6, tri, -35:5:30, true, "∂²f/∂x∂y")
    resize_to_layout!(fig)
    ε = 100sqrt(sum((norm.(to_mat.(H) .- f′′.(x, y))) .^ 2) / sum(norm.(to_mat.(H)) .^ 2))
    return fig, ε
end
_, Hg = generate_derivatives(tri, z)
fig, ε = plot_hessians(Hg, tri, f′′, x, y)
```

```julia-repl
julia> ε
42.085578794275605
```

```@raw html
<figure>
    <img src='../figures/hessian_data.png', alt'Cubic Hessians'><br>
</figure>
```

The error is certainly quite large, but remember that we are doing derivative _generation_ here rather than _estimation_. Judging from the figures themselves, the derivatives we have obtained are actually pretty good. 

Let's now see what happens if we only go up to quadratic terms.

```julia
_, Hg = generate_derivatives(tri, z, use_cubic_terms=false)
fig, ε = plot_hessians(Hg, tri, f′′, x, y)
```

```julia-repl
julia> ε
35.20873081559232
```

```@raw html
<figure>
    <img src='../figures/hessian_data_no_cubic.png', alt'Quadratic Hessians'><br>
</figure>
```

The error has actually decreased, and the figures do indeed look better. So, in this case, including cubic terms does not improve the results significantly (sometimes it does).

What if we used the iterative approach? 

```julia
_, Hg = generate_derivatives(tri, z, method=Iterative()) # the gradients will be generated first automatically
fig, ε = plot_hessians(Hg, tri, f′′, x, y)
```

```julia-repl
julia> ε
39.58481626576425
```

```@raw html
<figure>
    <img src='../figures/hessian_data_iterative.png', alt'Iterative Hessians'><br>
</figure>
```

The results are slightly worse, and varying `alpha` doesn't seem to do much.

# Generation Away from the Data Sites 

Now let's consider derivative generation away from the data sites. The function `differentiate` is used for this. We first construct our interpolant, ensuring we set `derivatives=true` so that we get the gradients at the data sites first, and then we `differentiate`.

```julia
itp = interpolate(tri, z; derivatives=true, method = Direct(), use_cubic_terms=false)
∂ = differentiate(itp, 1)
```

The second argument specifies the order of the resulting derivatives. Since we specify order 1, we will get gradients $(\partial_xf,\partial_yf)$.

Let's now define the grid for differentiating.

```julia
xg = LinRange(0, 1, 500)
yg = LinRange(0, 1, 500)
_x = vec([x for x in xg, _ in yg])
_y = vec([y for _ in xg, y in yg])
```

We can now evaluate `∂`. To approximate the function values at each point, we will use the `Sibson(1)` method, since this will incorporate the gradient information. I would really like to eventually get Hiyoshi's $C^2$ interpolant, as discussed in [Section 3.2.7.3 here](https://kluedo.ub.rptu.de/frontdoor/deliver/index/docId/2104/file/diss.bobach.natural.neighbor.20090615.pdf), as this will also give us $C^2$ continuity at the derivative sites and thus give smoother surfaces (noting some complexity issues discussed in Section 6.5 of the linked thesis), but I've just not found the time to comprehend how to best implement it yet / digest the spline notation (see issue [#1](https://github.com/DanielVandH/NaturalNeighbours.jl/issues/1) if you are interested on contributing to this). Note also that, just as with the interpolation methods, it is best to give vectors to `∂`. Lastly, since we are evaluating away from the data sites, remember that the Sibson coordinates are now incorporated into the weights of the associated weighted least squares problem (that you could disable if you for some reason wanted to with `use_sibson_weight=false`).

```julia
∇g = ∂(_x, _y; interpolant_method = Sibson(1))
```

Let's now plot our gradients. Note that there are some `Inf` values in the computed gradients, and these correspond to points evaluated outside of the convex hull of our data sites.

```julia
function rrmserr(z, ẑ)
    num = 0.0
    den = 0.0
    for (zᵢ, ẑᵢ) in zip(z, ẑ)
        if all(isfinite, (zᵢ..., ẑᵢ...))
            num += norm(zᵢ .- ẑᵢ)^2
            den += norm(ẑᵢ)^2
        end
    end
    # num /= length(ẑ)
    return 100sqrt(num / den)
end
function plot_f2(fig, x, y, vals, title, i, levels, show_3d=true, zlabel="z")
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = contourf!(ax, x, y, vals, colormap=:viridis, extendhigh=:auto, levels=levels)
    if show_3d
        ax = Axis3(fig[2, i], xlabel="x", ylabel="y", zlabel=zlabel, width=600, height=600, title=" ", titlealign=:left, azimuth=0.49)
        surface!(ax, x, y, vals, color=vals, colormap=:viridis, colorrange=extrema(levels))
    end
    return c
end
function plot_gradients(∇g, f′, xg, yg)
    ∇g = reshape(∇g, (length(xg), length(yg)))
    ∇g1 = first.(∇g)
    ∇g2 = last.(∇g)
    ∇f = [f′(x, y) for x in xg, y in yg]
    fig = Figure(fontsize=36, size=(2400, 600))
    plot_f2(fig, xg, yg, ∇g1, "(a): ∂f̂/∂x", 1, -3.5:0.5:3.0, true, "∂f̂/∂x")
    plot_f2(fig, xg, yg, ∇g2, "(b): ∂f̂/∂y", 3, -3.5:0.5:3.0, true, "∂f̂/∂y")
    plot_f2(fig, xg, yg, first.(∇f), "(c): ∂f/∂x", 2, -3.5:0.5:3.0, true, "∂f/∂x")
    plot_f2(fig, xg, yg, last.(∇f), "(d): ∂f/∂y", 4, -3.5:0.5:3.0, true, "∂f/∂y")
    plot_f2(fig, xg, yg, norm.(collect.(∇g) .- ∇f), "(e): Gradient error", 5, 0:0.1:0.5, true, "|∇ε|")
    resize_to_layout!(fig)
    ε = rrmserr(∇f, collect.(∇g))
    return fig, ε
end
fig, ε = plot_gradients(∇g, f′, xg, yg)
```

```julia-repl
julia> ε
13.185747607565729
```

```@raw html
<figure>
    <img src='../figures/gradient_surface.png', alt'Evaluated Gradient'><br>
</figure>
```

There are of course some strange artifacts near the convex hull, but the results are not terrible. Let's see what happens to the error if we instead use the other interpolant methods.

```julia
other_methods = [Sibson(), Laplace(), Nearest(), Triangle(; allow_cache = true)]
∇gs = [∂(_x, _y; interpolant_method=method) for method in other_methods]
∇f = [f′(x, y) for x in xg, y in yg]
εs = [rrmserr(∇f, collect.(∇g)) for ∇g in ∇gs]
```

```julia-repl
julia> hcat(other_methods, εs)
4×2 Matrix{Any}:
 Sibson{0}()    28.6753
 Laplace{0}()   25.499
 Nearest{0}()   69.5744
 Triangle{0}()  27.7737
```

Of course, the other methods are much worse.

Now let's go up to second order.

```julia
function plot_hessians(H, f′′, xg, yg)
    H = reshape(H, (length(xg), length(yg)))
    H₁₁ = getindex.(H, 1)
    H₁₂ = getindex.(H, 3)
    H₂₂ = getindex.(H, 2)
    Hf = [f′′(x, y) for x in xg, y in yg]
    fig = Figure(fontsize=36, size=(2400, 600))
    plot_f2(fig, xg, yg, H₁₁, "(a): ∂²f̂/∂x²", 1, -35:5:30, true, "∂²f̂/∂x²")
    plot_f2(fig, xg, yg, H₂₂, "(c): ∂²f̂/∂y²", 3, -35:5:30, true, "∂²f̂/∂y²")
    plot_f2(fig, xg, yg, H₁₂, "(e): ∂²f̂/∂x∂y", 5, -35:5:30, true, "∂²f̂/∂x∂y")
    plot_f2(fig, xg, yg, getindex.(Hf, 1), "(b): ∂²f/∂x²", 2, -35:5:30, true, "∂²f/∂x²")
    plot_f2(fig, xg, yg, getindex.(Hf, 4), "(d): ∂²f/∂y²", 4, -35:5:30, true, "∂²f/∂y²")
    plot_f2(fig, xg, yg, getindex.(Hf, 2), "(f): ∂²f/∂x∂y", 6, -35:5:30, true, "∂²f/∂x∂y")
    resize_to_layout!(fig)
    ε = rrmserr(Hf, to_mat.(H))
    return fig, ε
end
∂ = differentiate(itp, 2)
∇Hg = ∂(_x, _y; interpolant_method=Sibson(1), method = Iterative())
∇g = first.(∇Hg)
Hg = last.(∇Hg)
zlims!(figH.content[4], -25, 25)
fig∇, ε∇ = plot_gradients(∇g, f′, xg, yg)
figH, εH = plot_hessians(Hg, f′′, xg, yg)
```

```julia-repl
julia> ε∇
19.07546882353911

julia> εH
51.1267212244942
```

```@raw html
<figure>
    <img src='../figures/gradient_surface_2.png', alt'Evaluated Gradient'><br>
</figure>
```

```@raw html
<figure>
    <img src='../figures/hessian_surface.png', alt'Evaluated Hessian'><br>
</figure>
```

The gradients actually look better here, despite the greater error, especially around the convex hull. The Hessians are a bit problematic around the convex hull especially, but we are really asking a lot of the interpolant to get Hessians unfortunately. 

Let's see if the direct approach can give us any improvements (the default is `Iterative()` since we have derivative information in the interpolant).

```julia
∇Hg = ∂(_x, _y; interpolant_method=Sibson(1), method=Direct())
∇g = first.(∇Hg)
Hg = last.(∇Hg)
fig∇, ε∇ = plot_gradients(∇g, f′, xg, yg)
figH, εH = plot_hessians(Hg, f′′, xg, yg)
zlims!(figH.content[4], -25, 25)
```

```julia-repl
julia> ε∇
9.853286514069882

julia> εH
46.7610990050276
```

```@raw html
<figure>
    <img src='../figures/gradient_surface_2_direct.png', alt'Evaluated Gradient'><br>
</figure>
```

```@raw html
<figure>
    <img src='../figures/hessian_surface_direct.png', alt'Evaluated Hessian'><br>
</figure>
```

Indeed, both the gradients and the Hessians appear to have improved, with some difficulties at the convex hull. Perhaps a better way to measure the error is to only include points that are away fro the convex hull. The following function can do this for us:

```julia
function rrmserr(z, ẑ, tri, x, y)
    num = 0.0
    den = 0.0 
    points = get_point(tri)
    ch = get_convex_hull_vertices(tri)
    for (zᵢ, ẑᵢ, xᵢ, yᵢ) in zip(z, ẑ, x, y)
        q = (xᵢ, yᵢ)
        δ = DelaunayTriangulation.distance_to_polygon(q, points, ch)
        if δ > 1e-4 && all(isfinite, (zᵢ..., ẑᵢ...))
            num += norm(zᵢ .- ẑᵢ)^2
            den += norm(ẑᵢ)^2
        end
    end
    # num /= length(ẑ)
    return 100sqrt(num / den)
end
```

If we instead use this metric, then:

```julia-repl
julia> ε∇_nohull = rrmserr(f′.(_x, _y), ∇g, ∂, _x, _y)
7.479964687679311

julia> εH_nohull = rrmserr(f′′.(_x, _y), to_mat.(Hg), ∂, _x, _y)
38.884740966379056
```

The errors are smaller, though not by much.