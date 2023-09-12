using ..NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using ReferenceTests
using StableRNGs
using LinearAlgebra

f = (x, y) -> 0.75 * exp(-((9 * x - 2)^2 + (9 * y - 2)^2) / 4) + 0.75 * exp(-(9 * x + 1)^2 / 49 - (9 * y + 1) / 10) + 0.5 * exp(-((9 * x - 7)^2 + (9 * y - 3)^2) / 4) - 0.2 * exp(-(9 * x - 4)^2 - (9 * y - 7)^2)
f′ = (x, y) -> [(exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * x - 72)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * x) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * x) / 2 - 63 / 2)) / 2 - (3 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10) * ((162 * x) / 49 + 18 / 49)) / 4
    (exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * y - 126)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * y) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * y) / 2 - 27 / 2)) / 2 - (27 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)) / 40]
f′′ = (x, y) -> [(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/98-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49)^2)/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)^2)/5 (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5
    (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5 (243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/400+(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*y)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*y)/2-27/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*y-126)^2)/5]

function plot_f(fig, x, y, vals, title, i, show_3d=true, zlabel="z")
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = contourf!(ax, x, y, vals, color=vals, colormap=:viridis, extendhigh=:auto)
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

fig = Figure(fontsize=50)
plot_f(fig, x, y, z, "(a): f", 1, true, "z")
plot_f(fig, x, y, ∇₁, "(b): ∂f/∂x", 2, true, "∂f/∂x")
plot_f(fig, x, y, ∇₂, "(c): ∂f/∂y", 3, true, "∂f/∂y")
plot_f(fig, x, y, H₁₁, "(d): ∂²f/∂x²", 4, true, "∂²f/∂x²")
plot_f(fig, x, y, H₂₂, "(f): ∂²f/∂y²", 5, true, "∂²f/∂y²")
plot_f(fig, x, y, H₁₂, "(e): ∂²f/∂x∂y", 6, true, "∂²f/∂x∂y")
resize_to_layout!(fig)

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "differentiation_exact_surfaces.png") fig

# The data set 
rng = StableRNG(9199)
x = rand(rng, 500)
y = rand(rng, 500)
z = f.(x, y)
tri = triangulate([x'; y'])
vorn = voronoi(tri)

fig = Figure(fontsize=50, resolution=(1800, 600))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", width=600, height=600, title="(a): Data and triangulation", titlealign=:left)
scatter!(ax, x, y, color=:black, markersize=9)
triplot!(ax, tri, color=:black, linewidth=2)
voronoiplot!(ax, vorn, strokecolor=:blue, color=(:white, 0.0))
xlims!(ax, 0, 1)
ylims!(ax, 0, 1)

ax = Axis3(fig[1, 2], xlabel="x", ylabel="y", zlabel="z", width=600, height=600, azimuth=0.25, title="(b): Function values", titlealign=:left)
triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
surface!(ax, x, y, z)
scatter!(ax, x, y, z, color=:black, markersize=9)

colgap!(fig.layout, 1, 75)
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "example_data.png") fig

# Generating gradients at the data sites 
function plot_f2(fig, x, y, vals, title, i, tri, levels, show_3d=true, zlabel="z")
    triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = tricontourf!(ax, x, y, vals, color=vals, triangulation=triangles', colormap=:viridis, extendhigh=:auto, levels=levels)
    if show_3d
        ax = Axis3(fig[2, i], xlabel="x", ylabel="y", zlabel=zlabel, width=600, height=600, title=" ", titlealign=:left, azimuth=0.49)
        mesh!(ax, hcat(x, y, vals), triangles, color=vals, colormap=:viridis, colorrange=extrema(levels))
    end
    return c
end
function plot_gradients(∇g, tri, f′, x, y)
    ∇g1 = first.(∇g)
    ∇g2 = last.(∇g)
    fig = Figure(fontsize=50, resolution=(2400, 600))
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
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "gradient_data.png") fig
@test ε ≈ 10.2511800

∇gr, _ = generate_derivatives(tri, z; method=Direct())
fig, ε = plot_gradients(∇gr, tri, f′, x, y)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "joint_gradient_data.png") fig
@test ε ≈ 7.7177915977

# Hessians 
to_mat(H::NTuple{3,Float64}) = [H[1] H[3]; H[3] H[2]]
function plot_hessians(H, tri, f′′, x, y)
    H₁₁ = getindex.(H, 1)
    H₁₂ = getindex.(H, 3)
    H₂₂ = getindex.(H, 2)

    fig = Figure(fontsize=50, resolution=(2400, 600))
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
@test ε ≈ 42.085578794
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "hessian_data.png") fig

_, Hg = generate_derivatives(tri, z, use_cubic_terms=false)
fig, ε = plot_hessians(Hg, tri, f′′, x, y)
@test ε ≈ 35.20873081
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "hessian_data_no_cubic.png") fig

_, Hg = generate_derivatives(tri, z, method=Iterative()) # the gradients will be generated first automatically
fig, ε = plot_hessians(Hg, tri, f′′, x, y)
@test ε ≈ 39.584816
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "hessian_data_iterative.png") fig

# Away from the data sites 
itp = interpolate(tri, z; derivatives=true, method=Direct(), use_cubic_terms=false)
∂ = differentiate(itp, 1)

xg = LinRange(0, 1, 500)
yg = LinRange(0, 1, 500)
_x = vec([x for x in xg, _ in yg])
_y = vec([y for _ in xg, y in yg])

∇g = ∂(_x, _y; interpolant_method=Sibson(1))

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
    c = contourf!(ax, x, y, vals, color=vals, colormap=:viridis, extendhigh=:auto, levels=levels)
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
    fig = Figure(fontsize=50, resolution=(2400, 600))
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
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "gradient_surface.png") fig
@test ε ≈ 13.1857476

other_methods = [Sibson(), Laplace(), Nearest(), Triangle()]
∇gs = [∂(_x, _y; interpolant_method=method) for method in other_methods]
∇f = [f′(x, y) for x in xg, y in yg]
εs = [rrmserr(∇f, collect.(∇g)) for ∇g in ∇gs]
@test εs ≈ [28.6753, 25.499, 69.5744, 27.7737] rtol = 0.5

# Second order 
function plot_hessians(H, f′′, xg, yg)
    H = reshape(H, (length(xg), length(yg)))
    H₁₁ = getindex.(H, 1)
    H₁₂ = getindex.(H, 3)
    H₂₂ = getindex.(H, 2)
    Hf = [f′′(x, y) for x in xg, y in yg]
    fig = Figure(fontsize=50, resolution=(2400, 600))
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
∇Hg = ∂(_x, _y; interpolant_method=Sibson(1), method=Iterative())
∇g = first.(∇Hg)
Hg = last.(∇Hg)
fig∇, ε∇ = plot_gradients(∇g, f′, xg, yg)
figH, εH = plot_hessians(Hg, f′′, xg, yg)
zlims!(figH.content[4], -25, 25)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "hessian_surface.png") figH by = psnr_equality(18)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "gradient_surface_2.png") fig∇ by = psnr_equality(18)
@test ε∇ ≈ 19.07546882353911 rtol = 1e-1
@test εH ≈ 51.1267212244942 rtol = 1e-1

∇Hg = ∂(_x, _y; interpolant_method=Sibson(1), method=Direct())
∇g = first.(∇Hg)
Hg = last.(∇Hg)
fig∇, ε∇ = plot_gradients(∇g, f′, xg, yg)
figH, εH = plot_hessians(Hg, f′′, xg, yg)
zlims!(figH.content[4], -25, 25)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "hessian_surface_direct.png") figH
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "gradient_surface_2_direct.png") fig∇
@test εH ≈ 46.7610990050276
@test ε∇ ≈ 9.853286514069882

function rrmserr(z, ẑ, ∂, x, y)
    tri = ∂.interpolant.triangulation
    num = 0.0
    den = 0.0
    points = get_points(tri)
    ch = get_convex_hull_indices(tri)
    for (zᵢ, ẑᵢ, xᵢ, yᵢ) in zip(z, ẑ, x, y)
        q = (xᵢ, yᵢ)
        δ = DelaunayTriangulation.distance_to_polygon(q, points, ch)
        if δ > 1e-1 && all(isfinite, (zᵢ..., ẑᵢ...))
            num += norm(zᵢ .- ẑᵢ)^2
            den += norm(ẑᵢ)^2
        end
    end
    # num /= length(ẑ)
    return 100sqrt(num / den)
end

εH_nohull = rrmserr(f′′.(_x, _y), to_mat.(Hg), ∂, _x, _y)
ε∇_nohull = rrmserr(f′.(_x, _y), ∇g, ∂, _x, _y)
@test ε∇_nohull ≈ 7.479964687673911
@test εH_nohull ≈ 38.884740966379056