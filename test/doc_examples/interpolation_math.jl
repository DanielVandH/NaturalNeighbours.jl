using ..NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using ReferenceTests
using StableRNGs

include(normpath(@__DIR__, "../.", "helper_functions", "test_functions.jl"))

## Example of a tessellation
points = [
    (0.0, 0.0), (-1.0, 1.0), (-0.5, 1.0), (0.0, 1.0), (0.5, 1.0), (1.0, 1.0),
    (1.0, 0.8), (1.0, 0.0), (1.0, -0.5), (1.0, -1.0),
    (0.1, -1.0), (-0.8, -1.0), (-1.0, -1.0),
    (-1.0, -0.7), (-1.0, -0.1), (-1.0, 0.6),
    (-0.1, -0.8), (0.2, -0.8),
    (-0.6, -0.4), (0.9, 0.0), (-0.5, 0.5), (-0.4, 0.6), (-0.1, 0.8)
]
tri = triangulate(points)
vorn = voronoi(tri)
fig, ax, sc = voronoiplot(vorn, axis=(width=400, height=400), markercolor=:red)
resize_to_layout!(fig)
fig

@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "example_tessellation.png") fig

## Nearest neighbour interpolant 
f, f′, f′′ = test_4()
x1, y1 = point_set_1()
z1 = f.(x1, y1)
xg = LinRange(0, 1, 50)
yg = LinRange(0, 1, 50)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
ze = f.(x, y) |> x -> reshape(x, length(xg), length(yg))
itp = interpolate(x1, y1, z1)
vals = itp(x, y, method=:nearest) |> x -> reshape(x, length(xg), length(yg))
fig = Figure(fontsize=36, size=(1700, 600))
ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(a): $f^{\text{NEAR}}$", titlealign=:left)
surface!(ax, xg, yg, vals)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
ax = Axis3(fig[1, 2], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(b): $f$", titlealign=:left)
surface!(ax, xg, yg, ze)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "fnear_example.png") fig

## Laplace interpolation 
rng = StableRNG(872973)
pts = [(0.0, 8.0), (0.0, 0.0), (14.0, 0.0),
    (14.0, 8.0), (4.0, 4.0), (10.0, 6.0),
    (6.0, 2.0), (12.0, 4.0), (0.0, 4.0),
    (2.5, 5.0), (7.0, 3.3), (4.5, 5.2),
    (13.0, 0.5), (12.0, 6.0), (8.5, 3.5),
    (0.5, 6.0), (1.5, 6.0), (3.5, 6.0),
    (0.5, 2.0), (2.5, 2.0), (2.5, 2.5),
    (9.0, 2.0), (8.5, 6.0), (4.0, 2.0)]
tri = triangulate(pts, randomise=false, delete_ghosts=false, rng=rng)
vorn = voronoi(tri)
q = (5.0, 4.0)
tri2 = deepcopy(tri)
add_point!(tri2, q, rng=rng)
vorn2 = voronoi(tri2)
V = get_polygon(vorn2, DelaunayTriangulation.num_points(tri2))
AX2 = get_area(vorn2, DelaunayTriangulation.num_points(tri2))

fig, ax, sc = voronoiplot(vorn, axis=(width=400, height=400), markercolor=:red, markersize=7, color=:white)
xlims!(ax, 3, 9)
ylims!(ax, 1.5, 7)
Vcoords = [get_polygon_point(vorn2, i) for i in V]
poly!(ax, Vcoords, color=(:blue, 0.2), strokewidth=2, strokecolor=:blue)
scatter!(ax, [q], color=:magenta, markersize=14)
resize_to_layout!(fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "new_tile.png") fig

## Laplace interpolant 
f, f′, f′′ = test_4()
x1, y1 = point_set_1()
z1 = f.(x1, y1)
xg = LinRange(0, 1, 50)
yg = LinRange(0, 1, 50)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
ze = f.(x, y) |> x -> reshape(x, length(xg), length(yg))
itp = interpolate(x1, y1, z1)
vals = itp(x, y, method=Laplace()) |> x -> reshape(x, length(xg), length(yg))
fig = Figure(fontsize=36, size=(1700, 600))
ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(a): $f^{\text{LAP}}$", titlealign=:left)
surface!(ax, xg, yg, vals)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
ax = Axis3(fig[1, 2], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(b): $f$", titlealign=:left)
surface!(ax, xg, yg, ze)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "flap_example.png") fig

## Sibson interpolant 
f, f′, f′′ = test_4()
x1, y1 = point_set_1()
z1 = f.(x1, y1)
xg = LinRange(0, 1, 50)
yg = LinRange(0, 1, 50)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
ze = f.(x, y) |> x -> reshape(x, length(xg), length(yg))
itp = interpolate(x1, y1, z1)
vals = itp(x, y, method=Sibson()) |> x -> reshape(x, length(xg), length(yg))
fig = Figure(fontsize=36, size=(1700, 600))
ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(a): $f^{\text{SIB}0}$", titlealign=:left)
surface!(ax, xg, yg, vals)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
ax = Axis3(fig[1, 2], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(b): $f$", titlealign=:left)
surface!(ax, xg, yg, ze)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "fsib0_example.png") fig

## Sibson-1 interpolant 
f, f′, f′′ = test_4()
x1, y1 = point_set_1()
z1 = f.(x1, y1)
xg = LinRange(0, 1, 50)
yg = LinRange(0, 1, 50)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
ze = f.(x, y) |> x -> reshape(x, length(xg), length(yg))
itp = interpolate(x1, y1, z1, derivatives=true)
vals = itp(x, y, method=Sibson(1)) |> x -> reshape(x, length(xg), length(yg))
fig = Figure(fontsize=36, size=(1700, 600))
ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(a): $f^{\text{SIB}1}$", titlealign=:left)
surface!(ax, xg, yg, vals)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
ax = Axis3(fig[1, 2], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(b): $f$", titlealign=:left)
surface!(ax, xg, yg, ze)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "fsib1_example.png") fig

## Triangle interpolant 
f, f′, f′′ = test_4()
x1, y1 = point_set_1()
z1 = f.(x1, y1)
xg = LinRange(0, 1, 50)
yg = LinRange(0, 1, 50)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
ze = f.(x, y) |> x -> reshape(x, length(xg), length(yg))
itp = interpolate(x1, y1, z1, derivatives=true)
vals = itp(x, y, method=Triangle(; allow_cache = false)) |> x -> reshape(x, length(xg), length(yg))
fig = Figure(fontsize=36, size=(1700, 600))
ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(a): $f^{\text{TRI}}$", titlealign=:left)
surface!(ax, xg, yg, vals)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
ax = Axis3(fig[1, 2], xlabel=L"x", ylabel=L"y", zlabel=L"z", title=L"(b): $f$", titlealign=:left)
surface!(ax, xg, yg, ze)
scatter!(ax, x1, y1, z1, color=:red, strokewidth=2, markersize=4)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "ftri_example.png") fig
