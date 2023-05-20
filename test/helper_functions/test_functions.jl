# https://hdl.handle.net/10945/35052
function test_1()
    f = (x, y) -> 0.75 * exp(-((9 * x - 2)^2 + (9 * y - 2)^2) / 4) + 0.75 * exp(-(9 * x + 1)^2 / 49 - (9 * y + 1) / 10) + 0.5 * exp(-((9 * x - 7)^2 + (9 * y - 3)^2) / 4) - 0.2 * exp(-(9 * x - 4)^2 - (9 * y - 7)^2)
    f′ = (x, y) -> [(exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * x - 72)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * x) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * x) / 2 - 63 / 2)) / 2 - (3 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10) * ((162 * x) / 49 + 18 / 49)) / 4
        (exp(-(9 * x - 4)^2 - (9 * y - 7)^2) * (162 * y - 126)) / 5 - (3 * exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4) * ((81 * y) / 2 - 9)) / 4 - (exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4) * ((81 * y) / 2 - 27 / 2)) / 2 - (27 * exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)) / 40]
    f′′ = (x, y) -> [(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/98-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49)^2)/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)^2)/5 (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5
        (27*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10)*((162*x)/49+18/49))/40+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*x)/2-9)*((81*y)/2-9))/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*x)/2-63/2)*((81*y)/2-27/2))/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*x-72)*(162*y-126))/5 (243*exp(-(9 * y) / 10 - (9 * x + 1)^2 / 49 - 1 / 10))/400+(162*exp(-(9 * x - 4)^2 - (9 * y - 7)^2))/5-(243*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4))/8-(81*exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4))/4+(3*exp(-(9 * x - 2)^2 / 4 - (9 * y - 2)^2 / 4)*((81*y)/2-9)^2)/4+(exp(-(9 * x - 7)^2 / 4 - (9 * y - 3)^2 / 4)*((81*y)/2-27/2)^2)/2-(exp(-(9 * x - 4)^2 - (9 * y - 7)^2)*(162*y-126)^2)/5]
    return f, f′, f′′
end
function test_2()
    f = (x, y) -> (1 / 9) * (tanh(9 * y - 9 * x) + 1)
    f′ = (x, y) -> [tanh(9 * x - 9 * y)^2 - 1
        1 - tanh(9 * x - 9 * y)^2]
    f′′ = (x, y) -> [-2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9) 2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9)
        2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9) -2*tanh(9 * x - 9 * y)*(9*tanh(9 * x - 9 * y)^2-9)]
    return f, f′, f′′
end
function test_3()
    f = (x, y) -> (1.25 + cos(5.4 * y)) / (6 * (1 + (3 * x - 1)^2))
    f′ = (x, y) -> [-((108 * x - 36) * (cos((27 * y) / 5) + 5 / 4)) / (6 * (3 * x - 1)^2 + 6)^2
        -(27 * sin((27 * y) / 5)) / (5 * (6 * (3 * x - 1)^2 + 6))]
    f′′ = (x, y) -> [(2*(108*x-36)^2*(cos((27 * y) / 5)+5/4))/(6*(3*x-1)^2+6)^3-(108*(cos((27 * y) / 5)+5/4))/(6*(3*x-1)^2+6)^2 (27*sin((27 * y) / 5)*(108*x-36))/(5*(6*(3*x-1)^2+6)^2)
        (27*sin((27 * y) / 5)*(108*x-36))/(5*(6*(3*x-1)^2+6)^2) -(729 * cos((27 * y) / 5))/(25*(6*(3*x-1)^2+6))]
    return f, f′, f′′
end
function test_4()
    f = (x, y) -> (1 / 3) * exp(-(81 / 16) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
    f′ = (x, y) -> [-(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16) * ((81 * x) / 8 - 81 / 16)) / 3
        -(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16) * ((81 * y) / 8 - 81 / 16)) / 3]
    f′′ = (x, y) -> [(exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16))/8 (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)*((81*y)/8-81/16))/3
        (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*x)/8-81/16)*((81*y)/8-81/16))/3 (exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16)*((81*y)/8-81/16)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 16 - (81 * (y - 1 / 2)^2) / 16))/8]
    return f, f′, f′′
end
function test_5()
    f = (x, y) -> (1 / 3) * exp(-(81 / 4) * ((x - 1 / 2)^2 + (y - 1 / 2)^2))
    f′ = (x, y) -> [-(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4) * ((81 * x) / 2 - 81 / 4)) / 3
        -(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4) * ((81 * y) / 2 - 81 / 4)) / 3]
    f′′ = (x, y) -> [(exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4))/2 (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)*((81*y)/2-81/4))/3
        (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*x)/2-81/4)*((81*y)/2-81/4))/3 (exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4)*((81*y)/2-81/4)^2)/3-(27*exp(-(81 * (x - 1 / 2)^2) / 4 - (81 * (y - 1 / 2)^2) / 4))/2]
    return f, f′, f′′
end
function test_6()
    f = (x, y) -> (1 / 9) * (64 - 81 * ((x - 1 / 2)^2 + (y - 1 / 2)^2))^(1 / 2) - 1 / 2
    f′ = (x, y) -> [-(162 * x - 81) / (18 * (64 - 81 * (y - 1 / 2)^2 - 81 * (x - 1 / 2)^2)^(1 / 2))
        -(162 * y - 81) / (18 * (64 - 81 * (y - 1 / 2)^2 - 81 * (x - 1 / 2)^2)^(1 / 2))]
    f′′ = (x, y) -> [-(162 * x - 81)^2/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))-9/(64-81*(y-1/2)^2-81*(x-1/2)^2)^(1/2) -((162 * x - 81) * (162 * y - 81))/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))
        -((162 * x - 81) * (162 * y - 81))/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2)) -(162 * y - 81)^2/(36*(64-81*(y-1/2)^2-81*(x-1/2)^2)^(3/2))-9/(64-81*(y-1/2)^2-81*(x-1/2)^2)^(1/2)]
    return f, f′, f′′
end
function point_set_1()
    A = [0.022703 -0.031021
        0.021701 0.257692
        0.001903 0.494360
        0.039541 0.699342
        0.031583 0.910765
        0.132419 0.050133
        0.125444 0.259297
        0.076758 0.417112
        0.062649 0.655223
        0.095867 0.914652
        0.264560 0.029294
        0.208899 0.266878
        0.171473 0.480174
        0.190921 0.687880
        0.230463 0.904651
        0.366317 0.039695
        0.383239 0.238955
        0.346632 0.490299
        0.387316 0.644523
        0.379536 0.893803
        0.414977 -0.028462
        0.420001 0.226247
        0.485566 0.389142
        0.479258 0.632425
        0.397776 0.848971
        0.053989 0.158674
        0.017513 0.341401
        0.050968 0.578285
        0.048706 0.747019
        0.041878 0.996289
        0.109027 0.091855
        0.093454 0.338159
        0.145187 0.561556
        0.145273 0.752407
        0.069556 0.963242
        0.239164 0.060230
        0.276733 0.369604
        0.226678 0.594059
        0.186765 0.818558
        0.242622 0.980541
        0.385766 0.068448
        0.317909 0.312413
        0.377659 0.519930
        0.381292 0.820379
        0.280351 0.971172
        0.427768 0.156096
        0.466363 0.317509
        0.409203 0.508495
        0.481228 0.751101
        0.402732 0.997873
        0.584869 -0.027195
        0.606389 0.270927
        0.574131 0.425942
        0.599010 0.673378
        0.609697 0.924241
        0.661693 0.025596
        0.639647 0.200834
        0.700118 0.489070
        0.690895 0.669783
        0.671889 0.936610
        0.773694 0.028537
        0.741042 0.193658
        0.730603 0.471423
        0.821453 0.668505
        0.807664 0.847679
        0.842457 0.038050
        0.836692 0.208309
        0.847812 0.433563
        0.917570 0.630738
        0.927987 0.904231
        1.044982 -0.012090
        0.985788 0.269584
        1.012929 0.439605
        1.001985 0.694152
        1.041468 0.868208
        0.573008 0.127243
        0.501389 0.347773
        0.610695 0.608471
        0.538062 0.723524
        0.502619 1.030876
        0.642784 0.070783
        0.670396 0.325984
        0.633359 0.509632
        0.689564 0.775957
        0.683767 1.006451
        0.763533 0.102140
        0.825898 0.323577
        0.808661 0.609159
        0.729064 0.802281
        0.817095 1.051236
        0.868405 0.090205
        0.941846 0.331849
        0.859958 0.591014
        0.859633 0.814484
        0.851280 0.969603
        0.967063 0.133411
        0.967631 0.379528
        0.965704 0.504442
        1.035930 0.745992
        0.947151 0.980141]
    return A[:, 1], A[:, 2]
end
function point_set_2()
    A = [
        0.00 0.00
        0.00 1.00
        0.00 0.50
        0.50 1.00
        0.10 0.15
        0.15 0.30
        0.30 0.35
        0.10 0.75
        0.05 0.45
        1.00 0.00
        1.00 1.00
        0.50 0.00
        1.00 0.50
        0.20 0.10
        0.25 0.20
        0.60 0.25
        0.90 0.35
        0.80 0.40
        0.70 0.20
        0.95 0.90
        0.60 0.65
        0.65 0.70
        0.35 0.85
        0.60 0.85
        0.90 0.80
        0.85 0.25
        0.80 0.65
        0.75 0.85
        0.70 0.90
        0.70 0.65
        0.75 0.10
        0.75 0.35
        0.55 0.95
    ]
    return A[:, 1], A[:, 2]
end
function point_set_3()
    A = [
        0.1375 0.97500
        0.9125 0.98750
        0.7125 0.76250
        0.2250 0.83750
        0.0500 0.41250
        0.4750 0.63750
        0.0500 -0.05000
        0.4500 1.03750
        1.0875 0.55000
        0.5375 0.80000
        0.0375 0.75000
        0.1875 0.57500
        0.7125 0.55000
        0.8500 0.43750
        0.7000 0.31250
        0.2750 0.42500
        0.4500 0.28750
        0.8125 0.18750
        0.4500 -0.03750
        1.0000 0.26250
        0.5000 0.46250
        0.1875 0.26250
        0.5875 0.12500
        1.0500 -0.06125
        0.1000 0.11250
    ]
    return A[:, 1], A[:, 2]
end

function test_interpolant(itp, x, y, f)
    for method in (:sibson, :triangle, :nearest, :laplace, Sibson(1))
        for _ in 1:500
            vals = itp(x, y; parallel=false, method)
            vals2 = similar(vals)
            itp(vals2, x, y; parallel=false, method)
            vals3 = itp(x, y; parallel=true, method)
            vals4 = similar(vals3)
            itp(vals4, x, y; parallel=true, method)
            for i in eachindex(x, y)
                _x = x[i]
                _y = y[i]
                if method ≠ :nearest
                    _z = f isa Function ? f(_x, _y) : f[i]
                else
                    m = DT.jump_to_voronoi_polygon(itp.triangulation, (_x, _y))
                    _z = f isa Function ? f(get_point(itp.triangulation, m)...) : f[m]
                end
                @test all(val -> isapprox(val, _z, rtol=1e-1), (itp(_x, _y; method), vals[i], vals2[i], vals3[i], vals4[i]))
            end
        end
    end
end

function rrmse(z, ẑ)
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

function plot_2d_vals(fig, method, xp, yp, vals, errs, ε, i, ij, a, b, c, d, xg, yg, levels, x1, x2, title, show_scatter, x, y, z1, z2, z3)
    _i, _j = ij[i]
    ax = Axis(fig[_i, _j], xlabel=L"x", ylabel=L"y",
        title=title, titlealign=:left,
        width=400, height=400,
        xticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]), yticks=(0:0.25:1, [L"%$s" for s in 0:0.25:1]))
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    m = contourf!(ax, xg, yg, reshape(vals, length(xg), length(yg)), levels=levels, extendlow=:auto, extendhigh=:auto)
    if show_scatter
        color = xp == x1 ? :red : (xp == x2 ? :cyan : :magenta)
        scatter!(ax, xp, yp, color=color, markersize=7, strokecolor=:black)
    end
    i += 1
    return i, m
end

function plot_3d_vals(fig, method, xp, yp, vals, errs, ε, i, ij, a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, title, show_scatter, x, y, z1, z2, z3, zlabel)
    _i, _j = ij[i]
    ax = Axis3(fig[_i, _j], xlabel=L"x", ylabel=L"y", zlabel=zlabel,
        title=title, titlealign=:left,
        width=400, height=400,
        xticks=0:0.5:1, yticks=0:0.5:1, azimuth=0.9)
    xlims!(ax, a, b)
    ylims!(ax, c, d)
    zlims!(ax, zmin, zmax)
    tri = triangulate([x'; y'])
    triangles = [T[j] for T in each_solid_triangle(tri), j in 1:3]
    xyz = hcat(x, y, vals)
    m = mesh!(ax, xyz, triangles, color=vals, colorrange=(zmin, zmax))
    if show_scatter
        color = xp == x1 ? :red : (xp == x2 ? :cyan : :magenta)
        scatter!(ax, xp, yp, xp == x1 ? z1 : (xp == x2 ? z2 : z3), color=color, markersize=7, strokewidth=2)
    end
    i += 1
    return i, m
end

function update_dict_map(method, i, ε, dict_map, yval1, yval2, yval3, xtick, xticklabel, x1, x2, xp)
    if haskey(dict_map, method)
        i = dict_map[method]
    else
        if isempty(dict_map)
            dict_map[method] = i
        else
            i = maximum(values(dict_map)) + 1
            dict_map[method] = i
        end
    end
    if xp == x1
        push!(yval1, (ε, i))
    elseif xp == x2
        push!(yval2, (ε, i))
    else
        push!(yval3, (ε, i))
    end
    if method ∉ xticklabel
        push!(xtick, i)
        push!(xticklabel, string(method))
    end
    return i
end

function plot_scalar_summaries(fig, ii, j, method, ε, title, x1, x2, xps)
    xtick = Float64[]
    xticklabel = String[]
    yval1 = Tuple{Float64,Float64}[]
    yval2 = Tuple{Float64,Float64}[]
    yval3 = Tuple{Float64,Float64}[]
    dict_map = Dict()
    i = 1
    for (method, ε, xp) in zip(method, ε, xps)
        i = update_dict_map(method, i, ε, dict_map, yval1, yval2, yval3, xtick, xticklabel, x1, x2, xp)
    end
    unique!(xtick)
    unique!(xticklabel)
    ax = Axis(fig[ii, j], xlabel="RRMSE (%)", ylabel=" ", yticks=(xtick, xticklabel), title=title, titlealign=:left, width=600, height=400)
    scatter!(ax, yval1, color=:red, markersize=16, strokewidth=2)
    scatter!(ax, yval2, color=:cyan, markersize=16, strokewidth=2)
    scatter!(ax, yval3, color=:magenta, markersize=16, strokewidth=2)
    return fig
end

function complete_test_function_analysis(id)
    if id == 1
        f, f′, f′′ = test_1()
    elseif id == 2
        f, f′, f′′ = test_2()
    elseif id == 3
        f, f′, f′′ = test_3()
    elseif id == 4
        f, f′, f′′ = test_4()
    elseif id == 5
        f, f′, f′′ = test_5()
    elseif id == 6
        f, f′, f′′ = test_6()
    end
    (x1, y1), (x2, y2), (x3, y3) = point_set_1(), point_set_2(), point_set_3()
    xyt = [(x1, y1), (x2, y2), (x3, y3)]
    z1, z2, z3 = f.(x1, y1), f.(x2, y2), f.(x3, y3)
    zt = (z1, z2, z3)
    ∇1, ∇2, ∇3 = f′.(x1, y1), f′.(x2, y2), f′.(x3, y3)
    ∇t = (∇1, ∇2, ∇3)
    H1, H2, H3 = f′′.(x1, y1), f′′.(x2, y2), f′′.(x3, y3)
    Ht = (H1, H2, H3)
    xg = LinRange(0, 1, 50)
    yg = LinRange(0, 1, 50)
    x = vec([x for x in xg, _ in yg])
    y = vec([y for _ in xg, y in yg])
    ze, ∇e, He = f.(x, y), f′.(x, y), f′′.(x, y)
    interpolant_methods = (Sibson(), (Sibson(1), Direct()), (Sibson(1), Iterative()), Laplace(), Nearest(), Triangle())
    derivative_methods = (Direct(), Iterative())

    itp1 = interpolate(x1, y1, z1; derivatives=true)
    itp2 = interpolate(x2, y2, z2; derivatives=true)
    itp3 = interpolate(x3, y3, z3; derivatives=true)
    ∂11 = differentiate(itp1, 1)
    ∂12 = differentiate(itp1, 2)
    ∂21 = differentiate(itp2, 1)
    ∂22 = differentiate(itp2, 2)
    ∂31 = differentiate(itp3, 1)
    ∂32 = differentiate(itp3, 2)

    ## Get the errors
    # interpolant_results: Dict{Tuple{InterpolantMethod, PointSet}, {Values, Errors, RRMSE}}
    # gradient_results: Dict{Tuple{DifferentiatorMethod, InterpolantMethod, PointSet}, {Values, Errors, RRMSE}}
    # hessian_results: Dict{Tuple{DifferentiatorMethod, InterpolantMethod, PointSet}, {Values, Errors, RRMSE}}
    # data_site_gradient_results: Dict{Tuple{DifferentiatorMethod, PointSet}, {Values, Errors, RRMSE}}
    # data_site_hessian_results: Dict{Tuple{DifferentiatorMethod, PointSet}, {Values, Errors, RRMSE}}
    interpolant_results = OrderedDict{Any,Any}()
    first_order_gradient_results = OrderedDict{Any,Any}()
    second_order_gradient_results = OrderedDict{Any,Any}()
    second_order_hessian_results = OrderedDict{Any,Any}()
    gradient_data_site_generated_results = OrderedDict{Any,Any}()#
    first_order_gradient_data_site_evaluated_results = OrderedDict{Any,Any}()#
    second_order_gradient_data_site_evaluated_results = OrderedDict{Any,Any}()#
    hessian_data_site_generated_results = OrderedDict{Any,Any}()#
    second_order_hessian_data_site_evaluated_results = OrderedDict{Any,Any}()#

    # Interpolation errors
    for method in interpolant_methods
        for (i, itp) in pairs((itp1, itp2, itp3))
            if !(method isa Tuple)
                vals = itp(x, y; method=method)
            else
                _itp = interpolate(xyt[i][1], xyt[i][2], zt[i]; derivatives=true, method=method[2])
                vals = _itp(x, y; method=method[1])
            end
            all_errs = norm.(vals .- ze)
            err = rrmse(ze, vals)
            interpolant_results[(method, xyt[i])] = (vals, all_errs, err)
        end
    end

    # First order gradient errors 
    for method in derivative_methods
        for interpolant_method in interpolant_methods
            for (i, ∂) in pairs((∂11, ∂21, ∂31))
                if !(interpolant_method isa Tuple)
                    ∇ = ∂(x, y; method=method, interpolant_method=interpolant_method)
                else
                    _itp = interpolate(xyt[i][1], xyt[i][2], zt[i]; derivatives=true, method=interpolant_method[2])
                    _∂ = differentiate(_itp, 1)
                    ∇ = _∂(x, y; method=method, interpolant_method=interpolant_method[1])
                end
                all_errs = norm.(collect.(∇) .- ∇e)
                err = rrmse(∇e, ∇)
                first_order_gradient_results[(method, interpolant_method, xyt[i])] = (∇, all_errs, err)
            end
        end
    end

    # Second order gradient and Hessian errors
    for method in derivative_methods
        for interpolant_method in interpolant_methods
            for (i, ∂) in pairs((∂12, ∂22, ∂32))
                if !(interpolant_method isa Tuple)
                    ∇H = ∂(x, y; method=method, interpolant_method=interpolant_method)
                else
                    _itp = interpolate(xyt[i][1], xyt[i][2], zt[i]; derivatives=true, method=interpolant_method[2])
                    _∂ = differentiate(_itp, 2)
                    ∇H = _∂(x, y; method=method, interpolant_method=interpolant_method[1])
                end
                ∇ = first.(∇H)
                H = last.(∇H)
                all_∇_errs = norm.(collect.(∇) .- ∇e)
                err_∇ = rrmse(∇e, ∇)
                second_order_gradient_results[(method, interpolant_method, xyt[i])] = (∇, all_∇_errs, err_∇)
                all_H_errs = norm.(to_mat.(H) .- He)
                err_H = rrmse(He, to_mat.(H))
                second_order_hessian_results[(method, interpolant_method, xyt[i])] = (H, all_H_errs, err_H)
            end
        end
    end

    # Generated gradient and Hessian errors
    for (∇, H, (x, y), z) in zip(∇t, Ht, xyt, zt)
        for method in derivative_methods
            itp = interpolate(x, y, z; derivatives=true, method=method)
            ∇_generated = NNI.get_gradient(itp)
            H_generated = NNI.get_hessian(itp)
            all_∇_errs = norm.(collect.(∇_generated) .- ∇)
            err_∇ = rrmse(collect.(∇), collect.(∇_generated))
            gradient_data_site_generated_results[(method, (x, y))] = (∇_generated, all_∇_errs, err_∇)
            all_H_errs = norm.(to_mat.(H_generated) .- to_mat.(H))
            err_H = rrmse(to_mat.(H), to_mat.(H_generated))
            hessian_data_site_generated_results[(method, (x, y))] = (H_generated, all_H_errs, err_H)
        end
    end

    # First order gradient errors evaluated at the data sites 
    for method in derivative_methods
        for interpolant_method in interpolant_methods
            for (i, ∂) in pairs((∂11, ∂21, ∂31))
                if !(interpolant_method isa Tuple)
                    ∇ = ∂(xyt[i][1], xyt[i][2]; method=method, interpolant_method=interpolant_method)
                else
                    _itp = interpolate(xyt[i][1], xyt[i][2], zt[i]; derivatives=true, method=interpolant_method[2])
                    _∂ = differentiate(_itp, 1)
                    ∇ = _∂(xyt[i][1], xyt[i][2]; method=method, interpolant_method=interpolant_method[1])
                end
                all_errs = norm.(collect.(∇) .- ∇t[i])
                err = rrmse(∇t[i], ∇)
                first_order_gradient_data_site_evaluated_results[(method, interpolant_method, xyt[i])] = (∇, all_errs, err)
            end
        end
    end

    # Second order gradient and Hessian errors evaluated at the data sites 
    for method in derivative_methods
        for interpolant_method in interpolant_methods
            for (i, ∂) in pairs((∂12, ∂22, ∂32))
                if !(interpolant_method isa Tuple)
                    ∇H = ∂(xyt[i][1], xyt[i][2]; method=method, interpolant_method=interpolant_method)
                else
                    _itp = interpolate(xyt[i][1], xyt[i][2], zt[i]; derivatives=true, method=interpolant_method[2])
                    _∂ = differentiate(_itp, 2)
                    ∇H = _∂(xyt[i][1], xyt[i][2]; method=method, interpolant_method=interpolant_method[1])
                end
                ∇ = first.(∇H)
                H = last.(∇H)
                all_∇_errs = norm.(collect.(∇) .- ∇t[i])
                err_∇ = rrmse(∇t[i], ∇)
                second_order_gradient_data_site_evaluated_results[(method, interpolant_method, xyt[i])] = (∇, all_∇_errs, err_∇)
                all_H_errs = norm.(to_mat.(H) .- to_mat.(Ht[i]))
                err_H = rrmse(to_mat.(Ht[i]), to_mat.(H))
                second_order_hessian_data_site_evaluated_results[(method, interpolant_method, xyt[i])] = (H, all_H_errs, err_H)
            end
        end
    end

    ### Interpolant results 
    ij = [(1, 1), (2, 1), (3, 1),
        (1, 2), (2, 2), (3, 2),
        (1, 3), (2, 3), (3, 3),
        (1, 4), (2, 4), (3, 4),
        (1, 5), (2, 5), (3, 5),
        (1, 6), (2, 6), (3, 6)]
    a, b, c, d, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    all_errs = Float64[]
    for ((method, (xp, yp)), (vals, errs, ε)) in interpolant_results
        a = min(a, minimum(xp))
        b = max(b, maximum(xp))
        c = min(c, minimum(yp))
        d = max(d, maximum(yp))
        zmin = min(zmin, minimum(vals))
        zmax = max(zmax, maximum(vals))
        append!(all_errs, errs)
    end
    zmin = min(zmin, minimum(ze))
    zmax = max(zmax, maximum(ze))
    levels = LinRange(zmin, zmax, 20)
    εmin = quantile(all_errs, 0.02)
    εmax = quantile(all_errs, 0.98)
    εlevels = LinRange(εmin, εmax, 20)

    ## The interpolants
    fig = Figure(fontsize=36)
    i = 1
    for ((method, (xp, yp)), (vals, errs, ε)) in interpolant_results
        i, _ = plot_2d_vals(fig, (), xp, yp, vals, errs, ε, i, ij, a, b, c, d, xg, yg, levels, x1, x2, string(method), true, x, y, z1, z2, z3)
    end
    m = []
    for i in 1:6
        i, _m = plot_2d_vals(fig, (), (), (), ze, (), (), i, ((4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6)), a, b, c, d, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:4, 7], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "interpolation_results_2d", "test_function_$(id)_interpolation_results.png") fig

    ## 3D Interpolant plots 
    fig = Figure(fontsize=24)
    i = 1
    for ((method, (xp, yp)), (vals, errs, ε)) in interpolant_results
        i, _ = plot_3d_vals(fig, (), xp, yp, vals, errs, ε, i, ij, a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, string(method), true, x, y, z1, z2, z3, L"z")
    end
    m = []
    for i in 1:6
        i, _m = plot_3d_vals(fig, (), (), (), ze, (), (), i, ((4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6)), a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3, L"z")
        push!(m, _m)
    end
    Colorbar(fig[1:4, 7], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "interpolation_results_3d", "test_function_$(id)_3d_interpolation_results.png") fig

    ## The interpolation errors 
    fig = Figure(fontsize=36)
    i = 1
    m = []
    for ((method, (xp, yp)), (vals, errs, ε)) in interpolant_results
        i, _m = plot_2d_vals(fig, (), xp, yp, errs, errs, ε, i, ij, a, b, c, d, xg, yg, εlevels, x1, x2, string(method), true, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:3, 7], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "interpolation_errors", "test_function_$(id)_interpolation_error_results.png") fig

    ### Gradient results 
    ij = [(1, 1), (2, 1), (3, 1),
        (1, 2), (2, 2), (3, 2),
        (1, 3), (2, 3), (3, 3),
        (1, 4), (2, 4), (3, 4),
        (1, 5), (2, 5), (3, 5),
        (1, 6), (2, 6), (3, 6),
        (1, 7), (2, 7), (3, 7),
        (1, 8), (2, 8), (3, 8),
        (1, 9), (2, 9), (3, 9),
        (1, 10), (2, 10), (3, 10),
        (1, 11), (2, 11), (3, 11),
        (1, 12), (2, 12), (3, 12)]
    a, b, c, d, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    all_errs = Float64[]
    all_∇ = Float64[]
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_results
        a = min(a, minimum(xp))
        b = max(b, maximum(xp))
        c = min(c, minimum(yp))
        d = max(d, maximum(yp))
        nΔ = filter(x -> !any(isnan, x), vals)
        append!(all_∇, norm.(nΔ))
        append!(all_errs, filter(!isnan, errs))
    end
    zmin = 0.0
    zmax = quantile(all_∇, 0.99)
    levels = LinRange(zmin, zmax, 20)
    εmin = quantile(all_errs, 0.02)
    εmax = quantile(all_errs, 0.98)
    εlevels = LinRange(εmin, εmax, 20)

    ## The gradients
    fig = Figure(fontsize=36)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_results
        i, _ = plot_2d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
    end
    m = []
    for i in 1:12
        i, _m = plot_2d_vals(fig, (), (), (), norm.(∇e), (), (), i, [(4, i) for i in 1:12], a, b, c, d, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "first_order_gradient_results_2d", "test_function_$(id)_first_order_gradient_results.png") fig

    ## 3D Gradient plots
    fig = Figure(fontsize=24)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_results
        i, _ = plot_3d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3, L"|\nabla|")
    end
    m = []
    for i in 1:12
        i, _m = plot_3d_vals(fig, (), (), (), norm.(∇e), (), (), i, [(4, i) for i in 1:12], a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3, L"|\nabla|")
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "first_order_gradient_results_3d", "test_function_$(id)_3d_first_order_gradient_results.png") fig

    ## The gradient errors 
    fig = Figure(fontsize=36)
    i = 1
    m = []
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_results
        i, _m = plot_2d_vals(fig, (), xp, yp, errs, errs, ε, i, ij, a, b, c, d, xg, yg, εlevels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:3, 13], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "first_order_gradient_errors", "test_function_$(id)_first_order_gradient_error_results.png") fig

    ### Second order gradient results 
    a, b, c, d, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    all_errs = Float64[]
    all_∇ = Float64[]
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_results
        a = min(a, minimum(xp))
        b = max(b, maximum(xp))
        c = min(c, minimum(yp))
        d = max(d, maximum(yp))
        nΔ = filter(x -> !any(isnan, x), vals)
        append!(all_∇, norm.(nΔ))
        append!(all_errs, filter(!isnan, errs))
    end
    zmin = 0.0
    zmax = quantile(all_∇, 0.99)
    levels = LinRange(zmin, zmax, 20)
    εmin = quantile(all_errs, 0.02)
    εmax = quantile(all_errs, 0.98)
    εlevels = LinRange(εmin, εmax, 20)

    ## The gradients
    fig = Figure(fontsize=36)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_results
        i, _ = plot_2d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
    end
    m = []
    for i in 1:12
        i, _m = plot_2d_vals(fig, (), (), (), norm.(∇e), (), (), i, [(4, i) for i in 1:12], a, b, c, d, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_gradient_results_2d", "test_function_$(id)_second_order_gradient_results.png") fig

    ## 3D Gradient plots
    fig = Figure(fontsize=24)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_results
        i, _ = plot_3d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3, L"|\nabla|")
    end
    m = []
    for i in 1:12
        i, _m = plot_3d_vals(fig, (), (), (), norm.(∇e), (), (), i, [(4, i) for i in 1:12], a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3, L"|\nabla|")
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_gradient_results_3d", "test_function_$(id)_3d_second_order_gradient_results.png") fig

    ## The gradient errors 
    fig = Figure(fontsize=36)
    i = 1
    m = []
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_results
        i, _m = plot_2d_vals(fig, (), xp, yp, errs, errs, ε, i, ij, a, b, c, d, xg, yg, εlevels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:3, 13], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_gradient_errors", "test_function_$(id)_second_order_gradient_error_results.png") fig

    ### Hessian results
    a, b, c, d, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    all_errs = Float64[]
    all_H = Float64[]
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_results
        a = min(a, minimum(xp))
        b = max(b, maximum(xp))
        c = min(c, minimum(yp))
        d = max(d, maximum(yp))
        nΔ = filter(x -> !any(isnan, x), vals)
        append!(all_H, norm.(nΔ))
        append!(all_errs, filter(!isnan, errs))
    end
    zmin = 0.0
    zmax = quantile(all_H, 0.99)
    levels = LinRange(zmin, zmax, 20)
    εmin = quantile(all_errs, 0.02)
    εmax = quantile(all_errs, 0.98)
    εlevels = LinRange(εmin, εmax, 20)

    ## The Hessian
    fig = Figure(fontsize=36)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_results
        i, _ = plot_2d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
    end
    m = []
    for i in 1:12
        i, _m = plot_2d_vals(fig, (), (), (), norm.(He), (), (), i, [(4, i) for i in 1:12], a, b, c, d, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_hessian_results_2d", "test_function_$(id)_second_order_hessian_results.png") fig

    ## 3D Hessian plots
    fig = Figure(fontsize=24)
    i = 1
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_results
        i, _ = plot_3d_vals(fig, (), xp, yp, norm.(vals), errs, ε, i, ij, a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3, L"|H|")
    end
    m = []
    for i in 1:12
        i, _m = plot_3d_vals(fig, (), (),(), norm.(He), (),(), i, [(4, i) for i in 1:12], a, b, c, d, zmin, zmax, xg, yg, levels, x1, x2, "Exact", false, x, y, z1, z2, z3, L"|H|")
        push!(m, _m)
    end
    Colorbar(fig[1:4, 13], m[1])
    resize_to_layout!(fig)
    fig
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_hessian_results_3d", "test_function_$(id)_3d_second_order_hessian_results.png") fig

    ## The Hessian errors 
    fig = Figure(fontsize=36)
    i = 1
    m = []
    for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_results
        i, _m = plot_2d_vals(fig, (), xp, yp, errs, errs, ε, i, ij, a, b, c, d, xg, yg, εlevels, x1, x2, string(derivative_method) * string(interpolant_method), true, x, y, z1, z2, z3)
        push!(m, _m)
    end
    Colorbar(fig[1:3, 13], m[1])
    resize_to_layout!(fig)
    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "second_order_hessian_errors", "test_function_$(id)_second_order_hessian_error_results.png") fig

    ### Summarise errors
    ## Interpolation
    method_ε = [(method, ε, xp) for ((method, (xp, yp)), (vals, errs, ε)) in interpolant_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = Figure(fontsize=36)
    fig = plot_scalar_summaries(fig, 1, 1, method, ε, "Interpolation RRMSE", x1, x2, xp)

    ## First order gradinet
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 1, 2, method, ε, "First order gradient RRMSE", x1, x2, xp)

    ## Second order gradient errors 
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 1, 3, method, ε, "Second order gradient RRMSE", x1, x2, xp)

    ## Second order Hessian errors 
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 2, 1, method, ε, "Hessian RRMSE", x1, x2, xp)

    ## Generated gradients at the data sites 
    method_ε = [(method, ε, xp) for ((method, (xp, yp)), (vals, errs, ε)) in gradient_data_site_generated_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 2, 2, method, ε, "Generated gradient RRMSE", x1, x2, xp)

    ## Generated Hessians at the data sites 
    method_ε = [(method, ε, xp) for ((method, (xp, yp)), (vals, errs, ε)) in hessian_data_site_generated_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 2, 3, method, ε, "Generated Hessian RRMSE", x1, x2, xp)

    ## Evaluated first order gradients at the data sites 
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in first_order_gradient_data_site_evaluated_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 3, 1, method, ε, "Evaluated first order gradient RRMSE", x1, x2, xp)

    ## Evaluated second order gradients at the data sites
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_gradient_data_site_evaluated_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 3, 2, method, ε, "Evaluated second order gradient RRMSE", x1, x2, xp)

    ## Evaluated Hessians at the data sites
    method_ε = [((derivative_method, interpolant_method), ε, xp) for ((derivative_method, interpolant_method, (xp, yp)), (vals, errs, ε)) in second_order_hessian_data_site_evaluated_results]
    method, ε, xp = getindex.(method_ε, 1), getindex.(method_ε, 2), getindex.(method_ε, 3)
    fig = plot_scalar_summaries(fig, 3, 3, method, ε, "Evaluated Hessian RRMSE", x1, x2, xp)

    resize_to_layout!(fig)

    @test_reference normpath(@__DIR__, "..", "test_functions", "figures", "summary", "test_function_$(id)_summary.png") fig
end