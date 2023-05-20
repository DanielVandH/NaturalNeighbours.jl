using ..NaturalNeighbours
using CairoMakie
using DelaunayTriangulation
using DelimitedFiles
using Downloads
using StableRNGs
using StatsBase
using ReferenceTests

## Download and setup the data
data_url = "https://gist.githubusercontent.com/DanielVandH/13687b0918e45a416a5c93cd52c91449/raw/a8da6cdc94859fd66bcff85a2307f0f9cd57a18c/data.txt"
boundary_url = "https://gist.githubusercontent.com/DanielVandH/13687b0918e45a416a5c93cd52c91449/raw/a8da6cdc94859fd66bcff85a2307f0f9cd57a18c/boundary.txt"
data_dir = Downloads.download(data_url)
boundary_dir = Downloads.download(boundary_url)
data = readdlm(data_dir, skipstart=6)
data[:, 3] ./= 1000.0 # m -> km
boundary = readdlm(boundary_dir, skipstart=6)
good_elevation = findall(≥(0), @view data[:, 3])
data = @views data[good_elevation, :]
data_sites = [(data[i, 1], data[i, 2]) for i in axes(data, 1)]
elevation_data = @views data[:, 3]
boundary_points = [(boundary[i, 1], boundary[i, 2]) for i in axes(boundary, 1)]

## Setup the data for plotting 
## Need to get the triangulation for tricontourf, and we need to identify indices in data_sites for boundary_points 
function nearest_tuple(q, data)
    δ = Inf
    nearest_idx = 0
    qx, qy = getxy(q)
    for (i, p) in pairs(data)
        px, py = getxy(p)
        δ₁ = (qx - px)^2 + (qy - py)^2
        if δ₁ < δ
            δ = δ₁
            nearest_idx = i
        end
    end
    return nearest_idx
end
function update_boundary(boundary_points, data_sites)
    boundary_nodes = map(q -> nearest_tuple(q, data_sites), boundary_points)
    unique!(boundary_nodes)
    push!(boundary_nodes, boundary_nodes[begin])
    reverse!(boundary_nodes) # so that the boundary is traversed clockwise
    boundary_points = data_sites[boundary_nodes]
    return boundary_points, data_sites, boundary_nodes
end
boundary_points, data_sites, boundary_nodes = update_boundary(boundary_points, data_sites)
tri = triangulate(data_sites, boundary_nodes=boundary_nodes)
triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]

## Downsample the data to use for the interpolant 
rng = StableRNG(123)
desample_idx = sample(rng, axes(data, 1), 5000, replace=false)
ds_data = data[desample_idx, :]
ds_data_sites = [(ds_data[i, 1], ds_data[i, 2]) for i in axes(ds_data, 1)]
ds_elevation_data = @views ds_data[:, 3]
ds_boundary_points, ds_data_sites, ds_boundary_nodes = update_boundary(boundary_points, ds_data_sites)
reverse!(ds_boundary_nodes) # so that the boundary is traversed clockwise
ds_tri = triangulate(ds_data_sites, boundary_nodes=ds_boundary_nodes)
ds_triangles = [T[j] for T in each_solid_triangle(ds_tri), j in 1:3]

## Plot the data
colorrange = (0, 4)
levels = LinRange(colorrange..., 40)
fig = Figure(fontsize=24)
ax1 = Axis3(fig[1, 1], xlabel="Longitude", ylabel="Latitude", zlabel="Elevation (km)", width=600, height=400, azimuth=0.9, title="(a): Original height data (n = $(length(elevation_data)))", titlealign=:left)
mesh!(ax1, data, triangles, color=elevation_data, colorrange=colorrange)
ax2 = Axis(fig[1, 2], xlabel="Longitude", ylabel="Latitude", width=600, height=400, title="(b): Original height data (n = $(length(elevation_data)))", titlealign=:left)
tf = tricontourf!(ax2, data[:, 1], data[:, 2], elevation_data, color=elevation_data, triangulation=triangles', levels=levels)

ax3 = Axis3(fig[2, 1], xlabel="Longitude", ylabel="Latitude", zlabel="Elevation (km)", width=600, height=400, azimuth=0.9, title="(c): Downsampled height data (n = $(length(ds_elevation_data)))", titlealign=:left)
mesh!(ax3, ds_data, ds_triangles, color=ds_elevation_data, colorrange=colorrange)
ax4 = Axis(fig[2, 2], xlabel="Longitude", ylabel="Latitude", width=600, height=400, title="(d): Downsampled height data (n = $(length(ds_elevation_data)))", titlealign=:left)
tricontourf!(ax4, ds_data[:, 1], ds_data[:, 2], ds_elevation_data, color=ds_elevation_data, triangulation=ds_triangles', levels=levels)
Colorbar(fig[1:2, 3], tf)
resize_to_layout!(fig)
# save(normpath(@__DIR__, "..", "docs", "src", "figures", "swiss_heights.png"), fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "swiss_heights.png") fig

## Define the interpolant 
interpolant = interpolate(ds_data_sites, ds_elevation_data; derivatives=true)

## Evaluate the interpolant 
a, b, c, d = DelaunayTriangulation.polygon_bounds(ds_data_sites, ds_boundary_nodes)
nx = 250
ny = 250
xg = LinRange(a, b, nx)
yg = LinRange(c, d, ny)
x = [xg[i] for i in 1:nx, j in 1:ny] |> vec
y = [yg[j] for i in 1:nx, j in 1:ny] |> vec

sibson_vals = interpolant(x, y; method=Sibson(), parallel=true)
sibson_1_vals = interpolant(x, y; method=Sibson(1), parallel=true)
laplace_vals = interpolant(x, y; method=Laplace(), parallel=true)
triangle_vals = interpolant(x, y; method=Triangle(), parallel=true)
nearest_vals = interpolant(x, y; method=Nearest(), parallel=true)

## Plot the results 
query_tri = triangulate([x'; y'])
query_triangles = [T[j] for T in each_solid_triangle(query_tri), j in 1:3]
function plot_results!(fig, i1, j1, i2, j2, x, y, xg, yg, vals, title1, title2, query_triangles, query_tri, a, b, c, d, e, f, nx, ny)
    ax = Axis3(fig[i1, j1], xlabel="Longitude", ylabel="Latitude", zlabel="Elevation (km)", width=600, height=400, azimuth=0.9, title=title1, titlealign=:left)
    m = mesh!(ax, hcat(x, y, vals), query_triangles, color=vals, colorrange=colorrange)
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    zlims!(ax, e, f)
    ax = Axis(fig[i2, j2], xlabel="Longitude", ylabel="Latitude", width=600, height=400, title=title2, titlealign=:left)
    contourf!(ax, xg, yg, reshape(vals, (nx, ny)), color=vals, levels=levels)
    lines!(ax, [get_point(query_tri, i) for i in get_convex_hull_indices(query_tri)], color=:red, linewidth=4, linestyle=:dash)
    lines!(ax, ds_boundary_points, color=:white, linewidth=4)
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    return m
end
function plot_results(sibson_vals, sibson_1_vals, laplace_vals, triangle_vals, nearest_vals, query_triangles, interpolant, a, b, c, d, e, f, nx, ny, data, triangles, elevation_data)
    fig = Figure(fontsize=24)
    m1 = plot_results!(fig, 1, 1, 1, 2, x, y, xg, yg, sibson_vals, "(a): Sibson", "(b): Sibson", query_triangles, interpolant.triangulation, a, b, c, d, e, f, nx, ny)
    m2 = plot_results!(fig, 1, 3, 1, 4, x, y, xg, yg, sibson_1_vals, "(c): Sibson-1", "(d): Sibson-1", query_triangles, interpolant.triangulation, a, b, c, d, e, f, nx, ny)
    m3 = plot_results!(fig, 2, 1, 2, 2, x, y, xg, yg, laplace_vals, "(e): Laplace", "(f): Laplace", query_triangles, interpolant.triangulation, a, b, c, d, e, f, nx, ny)
    m4 = plot_results!(fig, 2, 3, 2, 4, x, y, xg, yg, triangle_vals, "(g): Triangle", "(h): Triangle", query_triangles, interpolant.triangulation, a, b, c, d, e, f, nx, ny)
    m5 = plot_results!(fig, 3, 1, 3, 2, x, y, xg, yg, nearest_vals, "(i): Nearest", "(j): Nearest", query_triangles, interpolant.triangulation, a, b, c, d, e, f, nx, ny)
    ax = Axis3(fig[3, 3], xlabel="Longitude", ylabel="Latitude", zlabel="Elevation (km)", width=600, height=400, azimuth=0.9, title="(k): Original height data", titlealign=:left)
    mesh!(ax, data, triangles, color=elevation_data, colorrange=(0, 4))
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    zlims!(ax, e, f)
    ax = Axis(fig[3, 4], xlabel="Longitude", ylabel="Latitude", width=600, height=400, title="(ℓ): Original height data", titlealign=:left)
    tricontourf!(ax, data[:, 1], data[:, 2], elevation_data, color=elevation_data, triangulation=triangles', levels=levels)
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    Colorbar(fig[1:3, 5], m1)
    resize_to_layout!(fig)
    return fig
end
e, f = 0.0, 4.5
fig = plot_results(sibson_vals, sibson_1_vals, laplace_vals, triangle_vals, nearest_vals, query_triangles, interpolant, a, b, c, d, e, f, nx, ny, data, triangles, elevation_data)
# save(normpath(@__DIR__, "..", "docs", "src", "figures", "swiss_heights_interpolated.png"), fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "swiss_heights_interpolated.png") fig

## Plot the results with project = false 
sibson_vals_p = interpolant(x, y; method=Sibson(), parallel=true, project=false)
sibson_1_vals_p = interpolant(x, y; method=Sibson(1), parallel=true, project=false)
laplace_vals_p = interpolant(x, y; method=Laplace(), parallel=true, project=false)
triangle_vals_p = interpolant(x, y; method=Triangle(), parallel=true, project=false)
nearest_vals_p = interpolant(x, y; method=Nearest(), parallel=true, project=false)

fig = plot_results(sibson_vals_p, sibson_1_vals_p, laplace_vals_p, triangle_vals_p, nearest_vals_p, query_triangles, interpolant, a, b, c, d, e, f, nx, ny, data, triangles, elevation_data)
# save(normpath(@__DIR__, "..", "docs", "src", "figures", "swiss_heights_interpolated_projected.png"), fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "swiss_heights_interpolated_projected.png") fig

## Plot the results, replacing all points outside of the boundary 
exterior_idx = identify_exterior_points(x, y, ds_data_sites, ds_boundary_nodes)
sibson_vals_p[exterior_idx] .= NaN
sibson_1_vals_p[exterior_idx] .= NaN
laplace_vals_p[exterior_idx] .= NaN
triangle_vals_p[exterior_idx] .= NaN
nearest_vals_p[exterior_idx] .= NaN

fig = plot_results(sibson_vals_p, sibson_1_vals_p, laplace_vals_p, triangle_vals_p, nearest_vals_p, query_triangles, interpolant, a, b, c, d, e, f, nx, ny, data, triangles, elevation_data)
# save(normpath(@__DIR__, "..", "docs", "src", "figures", "swiss_heights_interpolated_projected_boundary.png"), fig)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "swiss_heights_interpolated_projected_boundary.png") fig

