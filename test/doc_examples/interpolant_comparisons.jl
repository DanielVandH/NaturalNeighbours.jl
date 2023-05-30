using ..NaturalNeighbours
using CairoMakie
using ReferenceTests
using StableRNGs
using DelaunayTriangulation
using StaticArrays
using LinearAlgebra
using DataFrames
using StatsBase
using BenchmarkTools
using AlgebraOfGraphics
const NNI = NaturalNeighbours

## Some methods and constants 
const itp_methods = (
    Sibson(0),
    Triangle(),
    Nearest(),
    Laplace(),
    Sibson(1),
    Farin(1),
    Hiyoshi(2)
)
const diff_methods = (
    Direct(),
    Iterative()
)
const itp_aliases = (:Sibson0, :Triangle, :Nearest, :Laplace, :Sibson1, :Farin, :Hiyoshi)
const diff_aliases = (:Direct, :Iterative)
const itp_alias_map = Dict(itp_methods .=> itp_aliases)
const diff_alias_map = Dict(diff_methods .=> diff_aliases)
const colors = Dict(itp_aliases .=> [:red, :blue, :green, :orange, :purple, :black, :brown])
const linestyles = Dict(diff_aliases .=> [:solid, :dashdotdot])
const line_elements = [
    LineElement(color=color,
        linewidth=22,
        linestyle=:solid) for color in values(colors)
]
const style_elements = [
    LineElement(color=:black,
        linewidth=22,
        linestyle=linestyle) for linestyle in values(linestyles)
]
const azimuths = [0.3, 0.8, 0.3, 0.6, 0.6, 0.6, 0.45]
rng = StableRNG(123)
xg = LinRange(0, 1, 25)
yg = LinRange(0, 1, 25)
x = vec([x for x in xg, _ in yg])
y = vec([y for _ in xg, y in yg])
xg2 = LinRange(0, 1, 250)
yg2 = LinRange(0, 1, 250)
xq = vec([x for x in xg2, _ in yg2])
yq = vec([y for _ in xg2, y in yg2])
tol = 1e-2
tri = triangulate([x'; y']; rng=rng)
triq = triangulate([xq'; yq']; rng=rng)
exterior_idx = identify_exterior_points(xq, yq, get_points(tri), get_convex_hull_indices(tri); tol=tol)
interior_idx = filter(∉(exterior_idx), eachindex(xq, yq))

## The test functions 
const f = [
    (x, y) -> 0.75 * exp(-((9 * x - 2)^2 + (9 * y - 2)^2) / 4) + 0.75 * exp(-(9 * x + 1)^2 / 49 - (9 * y + 1) / 10) + 0.5 * exp(-((9 * x - 7)^2 + (9 * y - 3)^2) / 4) - 0.2 * exp(-(9 * x - 4)^2 - (9 * y - 7)^2)
    (x, y) -> (1 / 9) * (tanh(9 * y - 9 * x) + 1)
    (x, y) -> (1.25 + cos(5.4 * y)) / (6 * (1 + (3 * x - 1)^2))
    (x, y) -> (1 / 3) * exp(-(81 / 16) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
    (x, y) -> (1 / 3) * exp(-(81 / 4) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
    (x, y) -> (1 / 9) * (64 - 81 * ((x - 1 / 2)^2 + (y - 1 / 2)^2))^(1 / 2) - 1 / 2
    (x, y) -> sin(27 * x * y) - exp(-(x - y)^2 / 4) * cos(13 * x - 13 * y)
]
const ∇f = [
    (x, y) -> @SVector[(exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * x - 72)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * x) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * x) / 2 - 63 / 2)) / 2 - (3 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10) * ((162 * x) / 49 + 18 / 49)) / 4
        (exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * y - 126)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * y) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * y) / 2 - 27 / 2)) / 2 - (27 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)) / 40]
    (x, y) -> @SVector[tanh(9 * x - 9 * y)^2 - 1
        1 - tanh(9 * x - 9 * y)^2]
    (x, y) -> @SVector[-((108 * x - 36) * (cos((27 * y) / 5) + 5 / 4)) / (6 * (3 * x - 1)^2 + 6)^2
        -(27 * sin((27 * y) / 5)) / (5 * (6 * (3 * x - 1)^2 + 6))]
    (x, y) -> @SVector[-(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16) * ((81 * x) / 8 - 81 / 16)) / 3
        -(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16) * ((81 * y) / 8 - 81 / 16)) / 3]
    (x, y) -> @SVector[-(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4) * ((81 * x) / 2 - 81 / 4)) / 3
        -(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4) * ((81 * y) / 2 - 81 / 4)) / 3]
    (x, y) -> @SVector[-(162 * x - 81) / (18 * (64 - 81 * (y - 1 / 2)^2 - 81 * (x - 1 / 2)^2)^(1 / 2))
        -(162 * y - 81) / (18 * (64 - 81 * (y - 1 / 2)^2 - 81 * (x - 1 / 2)^2)^(1 / 2))]
    (x, y) -> @SVector[27 * y * cos(27 * x * y) + 13 * exp(-(x - y)^2 / 4) * sin(13 * x - 13 * y) + exp(-(x - y)^2 / 4) * cos(13 * x - 13 * y) * (x / 2 - y / 2)
        27 * x * cos(27 * x * y) - 13 * exp(-(x - y)^2 / 4) * sin(13 * x - 13 * y) - exp(-(x - y)^2 / 4) * cos(13 * x - 13 * y) * (x / 2 - y / 2)]
]
const Hf = [
    (x, y) -> @SMatrix[(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/98-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49)^2)/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)^2)/5 (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5
        (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5 (243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/400+(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*y)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*y)/2-27/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*y-126)^2)/5]
    (x, y) -> @SMatrix[-2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9) 2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9)
        2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9) -2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9)]
    (x, y) -> @SMatrix[(2*(108*x-36)^2*(cos((27 * y) / 5)+5/4))/(6*(3*x-1)^2+6)^3-(108*(cos((27 * y) / 5)+5/4))/(6*(3*x-1)^2+6)^2 (27*sin((27 * y) / 5)*(108*x-36))/(5*(6*(3*x-1)^2+6)^2)
        (27*sin((27 * y) / 5)*(108*x-36))/(5*(6*(3*x-1)^2+6)^2) -(729 * cos((27 * y) / 5))/(25*(6*(3*x-1)^2+6))]
    (x, y) -> @SMatrix[(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16))/8 (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)*((81*y)/8-81/16))/3
        (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)*((81*y)/8-81/16))/3 (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*y)/8-81/16)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16))/8]
    (x, y) -> @SMatrix[(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4))/2 (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)*((81*y)/2-81/4))/3
        (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)*((81*y)/2-81/4))/3 (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*y)/2-81/4)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4))/2]
    (x, y) -> @SMatrix[-(162 * x - 81)^2/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))-9/(64-81*(y-1/2)^2-81*(x-1/2)^2)^(1/2) -((162 * x - 81) * (162 * y - 81))/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))
        -((162 * x - 81) * (162 * y - 81))/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2)) -(162 * y - 81)^2/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))-9/(64-81*(y-1/2)^2-81*(x-1/2)^2)^(1/2)]
    (x, y) -> @SMatrix[(339*exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y))/2-729*y^2*sin(27 * x * y)-exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y)*(x/2-y/2)^2-26*exp(-(x - y)^2 / 4)*sin(13 * x - 13 * y)*(x/2-y/2) 27*cos(27 * x * y)-(339*exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y))/2-729*x*y*sin(27 * x * y)+exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y)*(x/2-y/2)^2+26*exp(-(x - y)^2 / 4)*sin(13 * x - 13 * y)*(x/2-y/2)
        27*cos(27 * x * y)-(339*exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y))/2-729*x*y*sin(27 * x * y)+exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y)*(x/2-y/2)^2+26*exp(-(x - y)^2 / 4)*sin(13 * x - 13 * y)*(x/2-y/2) (339*exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y))/2-729*x^2*sin(27 * x * y)-exp(-(x - y)^2 / 4)*cos(13 * x - 13 * y)*(x/2-y/2)^2-26*exp(-(x - y)^2 / 4)*sin(13 * x - 13 * y)*(x/2-y/2)]
]

## Functions for assessing qualities
function normal_to_triangle(p₁, p₂, p₃, z₁, z₂, z₃)
    x₁, y₁ = getxy(p₁)
    x₂, y₂ = getxy(p₂)
    x₃, y₃ = getxy(p₃)
    Δ = x₁ * y₂ - x₂ * y₁ - x₁ * y₃ + x₃ * y₁ + x₂ * y₃ - x₃ * y₂
    s₁ = (y₂ - y₃) / Δ
    s₂ = (y₃ - y₁) / Δ
    s₃ = (y₁ - y₂) / Δ
    s₄ = (x₃ - x₂) / Δ
    s₅ = (x₁ - x₃) / Δ
    s₆ = (x₂ - x₁) / Δ
    α = s₁ * z₁ + s₂ * z₂ + s₃ * z₃
    β = s₄ * z₁ + s₅ * z₂ + s₆ * z₃
    ∇norm = sqrt(1 + α^2 + β^2)
    ∇ = @SVector[-α, -β, 1.0]
    return ∇ / ∇norm
end
function normal_to_triangle(tri, z, i, j, k)
    p₁, p₂, p₃ = get_point(tri, i, j, k)
    z₁, z₂, z₃ = z[i], z[j], z[k]
    return normal_to_triangle(p₁, p₂, p₃, z₁, z₂, z₃)
end

function ∠(v₁, v₂)
    # acos is not reliable: https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf, https://scicomp.stackexchange.com/a/27694/42528
    a = norm(v₁)
    b = norm(v₂)
    c = norm(v₁ - v₂)
    if a < b
        a, b = b, a
    end
    μ = if b ≥ c
        c - (a - b)
    else
        b - (a - c)
    end
    num = ((a - b) + c) * μ
    den = (a + (b + c)) * ((a - c) + b)
    θ = 2atan(sqrt(num / den))
    return θ
end
function ∠(tri, z, i, j, k)
    # Angle between pᵢpⱼ and pᵢpₖ
    p₁, p₂, p₃ = get_point(tri, i, j, k)
    z₁, z₂, z₃ = z[i], z[j], z[k]
    px, py = getxy(p₁)
    qx, qy = getxy(p₂)
    rx, ry = getxy(p₃)
    v₁ = @SVector[qx - px, qy - py, z₂ - z₁]
    v₂ = @SVector[rx - px, ry - py, z₃ - z₁]
    return ∠(v₁, v₂)
end

function average_normal_vector(tri, z, i)
    # Using the mean-weighted-angle formula: https://doi.org/10.1007/s00371-004-0271-1
    n = @SVector[0.0, 0.0, 0.0]
    neighbouring_edges = get_adjacent2vertex(tri, i)
    for (j, k) in neighbouring_edges
        if !DelaunayTriangulation.is_ghost_triangle(i, j, k)
            ψ = ∠(tri, z, i, j, k)
            n = n + ψ * normal_to_triangle(tri, z, i, j, k)
        end
    end
    return n / norm(n)
end

function compare_normal_vectors(tri, z, i, ∇f::Function)
    # Maybe this is similar to https://doi.org/10.1007/978-3-319-40548-3_19?
    # The description is so vague.
    p = get_point(tri, i)
    x, y = getxy(p)
    n̄̂ = average_normal_vector(tri, z, i)
    nx, ny = ∇f(x, y)
    n = @SVector[-nx, -ny, 1.0]
    n̂ = n / norm(n)
    return rad2deg(∠(n̄̂, n̂))
end
function compare_normal_vectors(tri, z, ∇f::Function, interior_idx)
    return [compare_normal_vectors(tri, z, i, ∇f) for i in interior_idx]
end

function compare_quantities(ŷ, y, interior_idx)
    ε = 2norm.(ŷ .- y) ./ norm.(ŷ .+ y)
    return to_unit(ε[interior_idx])
end
function to_unit(μ)
    return max.(μ, sqrt(eps(Float64)))
end
to_mat(H) = @SMatrix[H[1] H[3]; H[3] H[2]]

## The analysis function 
function analysis_function!(df, tri, triq, x, y, xq, yq, fidx, itp_method, diff_method, interior_idx)
    g = f[fidx]
    ∇g = ∇f[fidx]
    Hg = Hf[fidx]
    z = g.(x, y)
    itp = interpolate(tri, z; derivatives=true, method=diff_method)
    ∂ = differentiate(itp, 2)

    ẑ = itp(xq, yq; method=itp_method)
    ∇̂Ĥ = ∂(xq, yq; method=diff_method, interpolant_method=itp_method)
    ∇̂ = SVector{2,Float64}.(first.(∇̂Ĥ))
    Ĥ = to_mat.(last.(∇̂Ĥ))

    z = g.(xq, yq)
    ∇ = ∇g.(xq, yq)
    H = Hg.(xq, yq)
    εz = compare_quantities(ẑ, z, interior_idx)
    ε∇ = compare_quantities(∇̂, ∇, interior_idx)
    εH = compare_quantities(Ĥ, H, interior_idx)
    εn = compare_normal_vectors(triq, ẑ, ∇g, interior_idx)

    _df = DataFrame(
        :z_exact => z[interior_idx],
        :z_approx => ẑ[interior_idx],
        :∇_exact => ∇[interior_idx],
        :∇_approx => ∇̂[interior_idx],
        :H_exact => H[interior_idx],
        :H_approx => Ĥ[interior_idx],
        :z_error => εz,
        :∇_error => ε∇,
        :H_error => εH,
        :n_error => εn,
        :itp_method => itp_alias_map[itp_method],
        :diff_method => diff_alias_map[diff_method],
        :f_idx => fidx
    )
    append!(df, _df)
    return df
end
function analysis_function(tri, triq, x, y, xq, yq, interior_idx)
    df = DataFrame(
        z_exact=Float64[],
        z_approx=Float64[],
        ∇_exact=SVector{2,Float64}[],
        ∇_approx=SVector{2,Float64}[],
        H_exact=SMatrix{2,2,Float64}[],
        H_approx=SMatrix{2,2,Float64}[],
        z_error=Float64[],
        ∇_error=Float64[],
        H_error=Float64[],
        n_error=Float64[],
        itp_method=Symbol[],
        diff_method=Symbol[],
        f_idx=Int[]
    )
    for fidx in eachindex(f, ∇f, Hf)
        for itp_method in itp_methods
            for diff_method in diff_methods
                analysis_function!(df, tri, triq, x, y, xq, yq, fidx, itp_method, diff_method, interior_idx)
            end
        end
    end
    return df
end

## Complete analysis 
df = analysis_function(tri, triq, x, y, xq, yq, interior_idx)
gdf = groupby(df, [:f_idx, :itp_method, :diff_method])

## Plot the results 
const alph = join('a':'z')
fig = Figure(fontsize=64)
z_ax = [Axis(fig[i, 1], xlabel=L"\varepsilon", ylabel=L"F(\varepsilon)",
    title=L"(%$(alph[i])1): $z$ $\varepsilon$ for $f_{%$i}", titlealign=:left,
    width=600, height=400, xscale=log10) for i in eachindex(f, ∇f, Hf)]
∇_ax = [Axis(fig[i, 2], xlabel=L"\varepsilon", ylabel=L"F(\varepsilon)",
    title=L"(%$(alph[i])2): $\nabla$ $\varepsilon$ for $f_{%$i}", titlealign=:left,
    width=600, height=400, xscale=log10) for i in eachindex(f, ∇f, Hf)]
H_ax = [Axis(fig[i, 3], xlabel=L"\varepsilon", ylabel=L"F(\varepsilon)",
    title=L"(%$(alph[i])3): $H$ $\varepsilon$ for $f_{%$i}", titlealign=:left,
    width=600, height=400, xscale=log10) for i in eachindex(f, ∇f, Hf)]
n_ax = [Axis(fig[i, 4], xlabel=L"\varepsilon", ylabel=L"F(\varepsilon)",
    title=L"(%$(alph[i])4): $n$ $\varepsilon$ for $f_{%$i}", titlealign=:left,
    width=600, height=400) for i in eachindex(f, ∇f, Hf)]
f_ax = [Axis3(fig[i, 5], xlabel=L"x", ylabel=L"y", zlabel=L"f_{%$i}(x, y)",
    title=L"(%$(alph[i])5): $f_{%$i}$'s surface", titlealign=:left,
    width=600, height=400, azimuth=azimuths[i]) for i in eachindex(f, ∇f, Hf)]
xℓ = [
    (1e-5, 1.0) (1e-3, 1.0) (1e-2, 1.0) (0.0, 5.0)
    (1e-5, 1.0) (1e-2, 1.0) (1e-1, 1.0) (0.0, 5.0)
    (1e-6, 1e-1) (1e-5, 1.0) (1e-2, 1.0) (0.0, 2.0)
    (1e-6, 1e-1) (1e-4, 1e-1) (1e-2, 1.0) (0.0, 1.0)
    (1e-5, 1e-1) (1e-3, 1.0) (1e-2, 1.0) (0.0, 2.0)
    (1e-8, 1e-1) (1e-5, 1e-1) (1e-2, 1e-1) (0.0, 0.5)
    (1e-2, 1.0) (1e-2, 1.0) (1e-1, 1.0) (0.0, 15.0)
]
for i in eachindex(f)
    xlims!(z_ax[i], xℓ[i, 1]...)
    xlims!(∇_ax[i], xℓ[i, 2]...)
    xlims!(H_ax[i], xℓ[i, 3]...)
    xlims!(n_ax[i], xℓ[i, 4]...)
end
for (f_idx, itp_alias, diff_alias) in keys(gdf)
    _df = gdf[(f_idx, itp_alias, diff_alias)]
    clr = colors[itp_alias]
    ls = linestyles[diff_alias]
    _z_ax = z_ax[f_idx]
    _∇_ax = ∇_ax[f_idx]
    _H_ax = H_ax[f_idx]
    _n_ax = n_ax[f_idx]
    z_error = _df.z_error
    ∇_error = _df.∇_error
    H_error = _df.H_error
    n_error = _df.n_error
    ecdfplot!(_z_ax, z_error, color=clr, linestyle=ls, linewidth=7)
    ecdfplot!(_∇_ax, ∇_error, color=clr, linestyle=ls, linewidth=7)
    ecdfplot!(_H_ax, H_error, color=clr, linestyle=ls, linewidth=7)
    if itp_alias ≠ :Nearest
        ecdfplot!(_n_ax, n_error, color=clr, linestyle=ls, linewidth=7)
    end
end
for f_idx in eachindex(f)
    g = f[f_idx]
    fz = [g(x, y) for x in xg2, y in yg2]
    _f_ax = f_ax[f_idx]
    surface!(_f_ax, xg2, yg2, fz)
end
[Legend(
    fig[i:(i+1), 6],
    [line_elements, style_elements],
    [string.(keys(colors)), string.(keys(linestyles))],
    ["Interpolant", "Differentiator"],
    titlesize=78,
    labelsize=78,
    patchsize=(100, 30)
) for i in (1, 3, 5)]
resize_to_layout!(fig)
fig
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "interpolant_comparison.png") fig

## Visual comparisons
considered_itp = eachindex(itp_methods)
considered_fidx = eachindex(f)
fig = Figure(fontsize=72, resolution=(4800, 4900))
ax = [
    Axis3(fig[i, j],
        xlabel=L"x",
        ylabel=L"y",
        zlabel=L"f(x, y)",
        title=L"(%$(alph[i])%$(j)): ($f_{%$i}$, %$(itp_aliases[j]))",
        titlealign=:left,
        width=600,
        height=600,
        azimuth=azimuths[i]
    )
    for i in considered_fidx, j in considered_itp
]
for (j, i) in enumerate(considered_fidx)
    for (ℓ, k) in enumerate(considered_itp)
        _gdf = gdf[(i, itp_aliases[k], diff_aliases[1])]
        _ax = ax[j, ℓ]
        _z = _gdf.z_approx
        surface!(_ax, xq[interior_idx], yq[interior_idx], _z)
        xlims!(_ax, 0, 1)
        ylims!(_ax, 0, 1)
        hidedecorations!(_ax)
    end
end
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "3d_visual_interpolant_comparison.png") fig

function plot_errors(considered_fidx, considered_itp, gdf, interior_idx, error_type, colorranges)
    fig = Figure(fontsize=72)
    ax = [
        Axis(fig[i, j],
            xlabel=L"x",
            ylabel=L"y",
            title=L"(%$(alph[i])%$(j)): ($f_{%$i}$, %$(itp_aliases[j]))",
            titlealign=:left,
            width=600,
            height=600,
        )
        for i in considered_fidx, j in considered_itp
    ]
    for (j, i) in enumerate(considered_fidx)
        for (ℓ, k) in enumerate(considered_itp)
            _gdf = gdf[(i, itp_aliases[k], diff_aliases[1])]
            _ax = ax[j, ℓ]
            ε = _gdf[!, error_type]
            heatmap!(_ax, xq[interior_idx], yq[interior_idx], ε, colorrange=colorranges[j])
            xlims!(_ax, 0, 1)
            ylims!(_ax, 0, 1)
            hidedecorations!(_ax)
        end
    end
    resize_to_layout!(fig)
    fig
end
z_colorranges = [(1e-4, 0.01), (1e-5, 0.1), (1e-2, 0.02), (1e-3, 0.005), (1e-2, 0.05), (1e-3, 0.005), (1e-2, 0.5)]
fig = plot_errors(considered_fidx, considered_itp, gdf, interior_idx, :z_error, z_colorranges)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "2d_visual_interpolant_comparison_z_error.png") fig

∇_colorranges = [(1e-2, 0.2), (1e-2, 1.0), (1e-2, 0.01), (1e-3, 0.01), (1e-2, 0.05), (1e-3, 0.01), (1e-2, 0.25)]
fig = plot_errors(considered_fidx, considered_itp, gdf, interior_idx, :∇_error, ∇_colorranges)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "2d_visual_interpolant_comparison_grad_error.png") fig

H_colorranges = [(1e-1, 0.5), (0.2, 0.8), (1e-1, 0.2), (1e-2, 0.2), (1e-1, 0.25), (1e-2, 0.1), (1e-1, 0.8)]
fig = plot_errors(considered_fidx, considered_itp, gdf, interior_idx, :H_error, H_colorranges)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "2d_visual_interpolant_comparison_hess_error.png") fig

n_colorranges = [(0, 5), (0, 5), (0, 2), (0, 1), (0, 2), (0, 2.5), (0, 15)]
fig = plot_errors(considered_fidx, considered_itp, gdf, interior_idx, :n_error, n_colorranges)
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "2d_visual_interpolant_comparison_normal_error.png") fig

## Random analysis 
function rrmse(y, ŷ) # interior_indices already applied
    num = 0.0
    den = 0.0
    for (yᵢ, ŷᵢ) in zip(y, ŷ)
        if all(isfinite, (yᵢ..., ŷᵢ...))
            num += norm(yᵢ .- ŷᵢ)^2
            den += norm(ŷᵢ)^2
        end
    end
    return 100sqrt(num / den)
end
function median_edge_length(tri)
    lengths = zeros(DelaunayTriangulation.num_solid_edges(tri))
    for (k, (i, j)) in (enumerate ∘ each_solid_edge)(tri)
        p, q = get_point(tri, i, j)
        px, py = getxy(p)
        qx, qy = getxy(q)
        ℓ = sqrt((qx - px)^2 + (qy - py)^2)
        lengths[k] = ℓ
    end
    return median(lengths)
end

function random_analysis_function(nsamples, triq, xq, yq, tol, rng)
    npoints = rand(rng, 50:2500)
    xs = [rand(rng, 50) for _ in 1:nsamples]
    ys = [rand(rng, 50) for _ in 1:nsamples]
    tris = [triangulate(tuple.(x, y); rng) for (x, y) in zip(xs, ys)]
    [refine!(tri; max_points=npoints) for tri in tris]
    xs = [first.(get_points(tri)) for tri in tris]
    ys = [last.(get_points(tri)) for tri in tris]
    exterior_idxs = [identify_exterior_points(xq, yq, get_points(tri), get_convex_hull_indices(tri); tol=tol) for tri in tris]
    interior_idxs = [filter(∉(exterior_idx), eachindex(xq, yq)) for exterior_idx in exterior_idxs]
    median_lengths = [median_edge_length(tri) for tri in tris]
    sortidx = sortperm(median_lengths)
    [permute!(obj, sortidx) for obj in (xs, ys, tris, exterior_idxs, interior_idxs, median_lengths)]
    dfs = Channel{DataFrame}(nsamples)
    Base.Threads.@threads for i in 1:nsamples
        tri = tris[i]
        x = xs[i]
        y = ys[i]
        interior_idx = interior_idxs[i]
        put!(dfs, analysis_function(tri, triq, x, y, xq, yq, interior_idx))
        println("Processed simulation $i.")
    end
    close(dfs)
    dfs = collect(dfs)
    df = DataFrame(
        f_idx=Int64[],
        itp_method=Symbol[],
        diff_method=Symbol[],
        z_rrmse=Float64[],
        ∇_rrmse=Float64[],
        H_rrmse=Float64[],
        n_error_median=Float64[],
        median_edge_length=Float64[]
    )
    for (i, _df) in enumerate(dfs)
        _gdf = groupby(_df, [:f_idx, :itp_method, :diff_method])
        _cgdf = combine(_gdf,
            [:z_exact, :z_approx] => ((z_exact, z_approx) -> rrmse(z_exact, z_approx)) => :z_rrmse,
            [:∇_exact, :∇_approx] => ((∇_exact, ∇_approx) -> rrmse(∇_exact, ∇_approx)) => :∇_rrmse,
            [:H_exact, :H_approx] => ((H_exact, H_approx) -> rrmse(H_exact, H_approx)) => :H_rrmse,
            :n_error => median => :n_error_median)
        _cgdf[!, :median_edge_length] .= median_lengths[i]
        append!(df, _cgdf)
    end
    _gdf = groupby(df, [:f_idx, :itp_method, :diff_method])
    return _gdf
end

nsamples = 50
rng = StableRNG(998881)
tol = 1e-1
random_results = random_analysis_function(nsamples, triq, xq, yq, tol, rng)

fig = Figure(fontsize=64)
z_ax = [Axis(fig[i, 1], xlabel=L"$ $Median edge length", ylabel=L"$z$ error",
    title=L"(%$(alph[i])1): $f_{%$i}", titlealign=:left,
    width=600, height=400, yscale=log10) for i in eachindex(f, ∇f, Hf)]
∇_ax = [Axis(fig[i, 2], xlabel=L"$ $Median edge length", ylabel=L"$\nabla$ error",
    title=L"(%$(alph[i])2): $f_{%$i}", titlealign=:left,
    width=600, height=400, yscale=log10) for i in eachindex(f, ∇f, Hf)]
H_ax = [Axis(fig[i, 3], xlabel=L"$ $Median edge length", ylabel=L"$H$ error",
    title=L"(%$(alph[i])3): $f_{%$i}", titlealign=:left,
    width=600, height=400, yscale=log10) for i in eachindex(f, ∇f, Hf)]
n_ax = [Axis(fig[i, 4], xlabel=L"$ $Median edge length", ylabel=L"$n$ error",
    title=L"(%$(alph[i])4): $f_{%$i}", titlealign=:left,
    width=600, height=400, yscale=log10) for i in eachindex(f, ∇f, Hf)]
for (f_idx, itp_alias, diff_alias) in keys(random_results)
    _df = random_results[(f_idx, itp_alias, diff_alias)]
    _df = filter(:itp_method => !=(:Nearest), _df)
    clr = colors[itp_alias]
    ls = linestyles[diff_alias]
    _z_ax = z_ax[f_idx]
    _∇_ax = ∇_ax[f_idx]
    _H_ax = H_ax[f_idx]
    _n_ax = n_ax[f_idx]
    x = _df.median_edge_length
    z_error = _df.z_rrmse
    ∇_error = _df.∇_rrmse
    H_error = _df.H_rrmse
    n_error = _df.n_error_median
    lines!(_z_ax, x, z_error, color=clr, linestyle=ls, linewidth=7)
    lines!(_∇_ax, x, ∇_error, color=clr, linestyle=ls, linewidth=7)
    lines!(_H_ax, x, H_error, color=clr, linestyle=ls, linewidth=7)
    lines!(_n_ax, x, n_error, color=clr, linestyle=ls, linewidth=7)
end
[Legend(
    fig[i:(i+1), 6],
    [line_elements, style_elements],
    [string.(keys(colors)), string.(keys(linestyles))],
    ["Interpolant", "Differentiator"],
    titlesize=78,
    labelsize=78,
    patchsize=(100, 30)
) for i in (1, 3, 5)]
resize_to_layout!(fig)
fig
@test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "median_edge_length_comparisons.png") fig by = psnr_equality(10)

## Computation times 
#=
if get(ENV, "CI", "false") == "false"
    function circular_example(m) # extra points are added outside of the circular barrier for derivative generation
        pts = [(cos(θ) + 1e-6randn(), sin(θ) + 1e-6randn()) for θ = LinRange(0, 2π, (m + 1))][1:end-1] # avoid cocircular points
        extra_pts = NTuple{2,Float64}[]
        while length(extra_pts) < 50
            p = (5randn(), 5randn())
            if norm(p) > 1.01
                push!(extra_pts, p)
            end
        end
        append!(pts, extra_pts)
        tri = triangulate(pts)
        return tri
    end

    function running_time_analysis(itp_method, m_range, g)
        running_times = zeros(length(m_range))
        for (i, m) in enumerate(m_range)
            tri = circular_example(m)
            z = [g(x, y) for (x, y) in each_point(tri)]
            itp = interpolate(tri, z; derivatives=true)
            b = @benchmark $itp($0.0, $0.0; method=$itp_method)
            running_times[i] = minimum(b.times) / 1e6 # ms
        end
        return DataFrame(
            running_times=running_times,
            method=itp_alias_map[itp_method],
            m=m_range
        )
    end
    function running_time_analysis(m_range, g)
        df = DataFrame(
            running_times=Float64[],
            method=Symbol[],
            m=Int[]
        )
        for itp_method in itp_methods
            _running_times = running_time_analysis(itp_method, m_range, g)
            append!(df, _running_times)
        end
        return df
    end

    m_range = 3:20
    g = f[end]
    running_times = running_time_analysis(m_range, g)

    fig = data(running_times) *
          mapping(:m, :running_times) *
          mapping(color=:method) *
          visual(Scatter, markersize=14) |> plt ->
        draw(plt; axis=(width=600, height=400, yscale=log10, xlabel=L"$ $Number of natural neighbours", ylabel=L"$t$ (ms)"))
    vlines!(fig.figure[1, 1], [6], linewidth=3, linestyle=:dash, color=:black)
    @test_reference normpath(@__DIR__, "../..", "docs", "src", "figures", "method_benchmarks.png") fig by = psnr_equality(10)
end
=#