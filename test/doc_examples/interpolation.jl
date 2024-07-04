using ..NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using ReferenceTests
using StableRNGs

include(normpath(@__DIR__, "../.", "helper_functions", "test_functions.jl"))

f = (x, y) -> (1 / 9) * (tanh(9 * y - 9 * x) + 1)

rng = StableRNG(123)
x = rand(rng, 100)
y = rand(rng, 100)
z = f.(x, y)
itp = interpolate(x, y, z; derivatives=true)

xg = LinRange(0, 1, 100)
yg = LinRange(0, 1, 100)
_x = vec([x for x in xg, _ in yg])
_y = vec([y for _ in xg, y in yg])
sib_vals = itp(_x, _y)
sib1_vals = itp(_x, _y; method=Sibson(1))

function plot_itp(fig, x, y, vals, title, i, show_data_sites, itp, xd=nothing, yd=nothing, show_3d=true, levels=-0.1:0.05:0.3)
    ax = Axis(fig[1, i], xlabel="x", ylabel="y", width=600, height=600, title=title, titlealign=:left)
    c = contourf!(ax, x, y, vals, colormap=:viridis, levels=levels, extendhigh=:auto)
    show_data_sites && scatter!(ax, xd, yd, color=:red, markersize=9)
    tri = itp.triangulation
    ch_idx = get_convex_hull_vertices(tri)
    lines!(ax, [get_point(tri, i) for i in ch_idx], color=:white, linewidth=4)
    if show_3d
        ax = Axis3(fig[2, i], xlabel="x", ylabel="y", zlabel=L"z", width=600, height=600, title=" ", titlealign=:left, azimuth=0.49)
        surface!(ax, x, y, vals, color=vals, colormap=:viridis, colorrange=extrema(levels))
        zlims!(ax, 0, 0.25)
    end
    return c
end
fig = Figure(fontsize=36)
plot_itp(fig, _x, _y, sib_vals, "(a): Sibson", 1, false, itp, x, y)
plot_itp(fig, _x, _y, sib1_vals, "(b): Sibson-1", 2, false, itp, x, y)
plot_itp(fig, _x, _y, f.(_x, _y), "(c): Exact", 3, true, itp, x, y)
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "sibson_vs_sibson1.png") fig

sib_vals = itp(_x, _y, project=false)
sib1_vals = itp(_x, _y; method=Sibson(1), project=false)
fig = Figure(fontsize=36)
plot_itp(fig, _x, _y, sib_vals, "(a): Sibson", 1, false, itp, x, y)
plot_itp(fig, _x, _y, sib1_vals, "(b): Sibson-1", 2, false, itp, x, y)
plot_itp(fig, _x, _y, f.(_x, _y), "(c): Exact", 3, true, itp, x, y)
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "sibson_vs_sibson1_project_false.png") fig

sib_vals = itp(_x, _y)
sib1_vals = itp(_x, _y; method=Sibson(1))
sib_errs = @. 100abs(sib_vals - f(_x, _y))
sib1_errs = @. 100abs(sib1_vals - f(_x, _y))

fig = Figure(fontsize=36)
plot_itp(fig, _x, _y, sib_errs, "(a): Sibson", 1, true, itp, x, y, false, 0:0.5:3)
c = plot_itp(fig, _x, _y, sib1_errs, "(b): Sibson-1", 2, true, itp, x, y, false, 0:0.5:3)
Colorbar(fig[1, 3], c)
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "sibson_vs_sibson1_errors.png") fig by=psnr_equality(20)

esib0 = 100sqrt(sum((sib_vals .- f.(_x, _y)).^2) / sum(sib_vals.^2))
esib1 = 100sqrt(sum((sib1_vals .- f.(_x, _y)).^2) / sum(sib_vals.^2))

@test esib1 < esib0