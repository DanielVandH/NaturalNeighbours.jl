function objective_function_gradient_direct(θ, p)
    tri, S, z, xᵢ, yᵢ, zᵢ, λ = p
    ℓ = 0.0
    β₁, β₂ = θ
    for (λ, s) in zip(λ, S)
        zₛ = z[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        wᵢ = λ / ((xₛ - xᵢ)^2 + (yₛ - yᵢ)^2)
        εᵢ² = wᵢ * (zₛ - zᵢ - β₁ * (xₛ - xᵢ) - β₂ * (yₛ - yᵢ))^2
        ℓ += εᵢ²
    end
    return ℓ
end
function estimate_gradient_direct(tri, r, z; use_sibson_weight=false, rng=Random.default_rng())
    if r isa Integer
        S = DT.iterated_neighbourhood(tri, r, 1)
        zᵢ = z[r]
        xᵢ, yᵢ = getxy(get_point(tri, r))
        S = collect(S)
        λ = similar(S)
        fill!(λ, 1)
    else
        nc = NNI.compute_natural_coordinates(Sibson(), tri, r, rng=rng)
        S = NNI.get_indices(nc)
        λ = NNI.get_coordinates(nc)
        zᵢ = 0.0
        for (a, b) in zip(λ, S)
            zₖ = z[b]
            zᵢ += a * zₖ
        end
        xᵢ, yᵢ = getxy(r)
        if !use_sibson_weight
            fill!(λ, 1.0)
        end
    end

    r = (xᵢ, yᵢ)
    X1 = [getx(get_point(tri, s)) - getx(r) for s in S]
    X2 = [gety(get_point(tri, s)) - gety(r) for s in S]
    X = hcat(X1, X2)
    W = [sqrt(λ) / norm(r .- get_point(tri, s)) for (λ, s) in zip(λ, S)]
    W .= W .^ 2
    W = Diagonal(W)
    Z = [z[s] - zᵢ for s in S]
    ∇ = (X' * W * X) \ (X' * W * Z)

    prob = OptimizationProblem(
        objective_function_gradient_direct,
        copy(∇),
        (tri, S, z, xᵢ, yᵢ, zᵢ, λ)
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD)
    return sol.u, ∇
end

function objective_function_gradient_hessian_cubic_direct(θ, p)
    tri, S, z, xᵢ, yᵢ, zᵢ = p
    ℓ = 0.0
    ∇₁, ∇₂, H₁₁, H₂₂, H₁₂, Q₁, Q₂, Q₃, Q₄ = θ
    for s in S
        zₛ = z[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        wᵢ = 1 / ((xₛ - xᵢ)^2 + (yₛ - yᵢ)^2)
        term₁ = ∇₁ * (xₛ - xᵢ) + ∇₂ * (yₛ - yᵢ)
        term₂ = H₁₁ * (xₛ - xᵢ)^2 / 2 + H₂₂ * (yₛ - yᵢ)^2 / 2 + H₁₂ * (xₛ - xᵢ) * (yₛ - yᵢ)
        term₃ = Q₁ * (xₛ - xᵢ)^3 / 6 + Q₂ * (xₛ - xᵢ)^2 * (yₛ - yᵢ) / 2 + Q₃ * (xₛ - xᵢ) * (yₛ - yᵢ)^2 / 2 + Q₄ * (yₛ - yᵢ)^3 / 6
        taylor = zᵢ + term₁ + term₂ + term₃
        εᵢ² = wᵢ * (zₛ - taylor)^2
        ℓ += εᵢ²
    end
    return ℓ
end
function estimate_gradient_hessian_cubic_direct(tri, r, z; rng=Random.default_rng())
    if r isa Integer
        S = DT.iterated_neighbourhood(tri, r, 3)
        zᵢ = z[r]
        xᵢ, yᵢ = getxy(get_point(tri, r))
    else
        nc = NNI.compute_natural_coordinates(Sibson(), tri, r, rng=rng)
        S = NNI.get_indices(nc)
        λ = NNI.get_coordinates(nc)
        zᵢ = 0.0
        for (a, b) in zip(λ, S)
            zₖ = z[b]
            zᵢ += a * zₖ
        end
        xᵢ, yᵢ = getxy(r)
        _S = [get_neighbours(tri, i) for i in nc.indices]
        _S1 = copy(nc.indices)
        push!(_S1, reduce(union, _S)...)
        filter!(!DT.is_ghost_vertex, _S1)
        unique!(_S1)
        _S2 = [get_neighbours(tri, i) for i in _S1]
        push!(_S1, reduce(union, _S2)...)
        filter!(!DT.is_ghost_vertex, _S1)
        unique!(_S1)
        S = _S1
    end

    r = (xᵢ, yᵢ)
    x = [getx(get_point(tri, s)) - getx(r) for s in S]
    y = [gety(get_point(tri, s)) - gety(r) for s in S]
    X = @. [x y x^2 / 2 y^2 / 2 x * y x^3 / 6 y^3 / 6 x^2 * y / 2 x * y^2 / 2]
    W = [1 / norm(r .- get_point(tri, s)) for s in S]
    W .= W .^ 2
    W = Diagonal(W)
    Z = [z[s] - zᵢ for s in S]
    ∇ℋ = (X' * W * X) \ (X' * W * Z)

    prob = OptimizationProblem(
        objective_function_gradient_hessian_cubic_direct,
        copy(∇ℋ),
        (tri, S, z, xᵢ, yᵢ, zᵢ)
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD)

    return (sol.u[1:2], sol.u[3:5]), (∇ℋ[1:2], ∇ℋ[3:5])
end

function objective_function_gradient_hessian_quadratic_direct(θ, p)
    tri, S, z, xᵢ, yᵢ, zᵢ = p
    ℓ = 0.0
    ∇₁, ∇₂, H₁₁, H₂₂, H₁₂ = θ
    for s in S
        zₛ = z[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        wᵢ = 1 / ((xₛ - xᵢ)^2 + (yₛ - yᵢ)^2)
        term₁ = ∇₁ * (xₛ - xᵢ) + ∇₂ * (yₛ - yᵢ)
        term₂ = H₁₁ * (xₛ - xᵢ)^2 / 2 + H₂₂ * (yₛ - yᵢ)^2 / 2 + H₁₂ * (xₛ - xᵢ) * (yₛ - yᵢ)
        taylor = zᵢ + term₁ + term₂
        εᵢ² = wᵢ * (zₛ - taylor)^2
        ℓ += εᵢ²
    end
    return ℓ
end
function estimate_gradient_hessian_quadratic_direct(tri, r, z; rng=Random.default_rng())
    if r isa Integer
        S = DT.iterated_neighbourhood(tri, r, 2)
        zᵢ = z[r]
        xᵢ, yᵢ = getxy(get_point(tri, r))
    else
        nc = NNI.compute_natural_coordinates(Sibson(), tri, r, rng=rng)
        S = NNI.get_indices(nc)
        λ = NNI.get_coordinates(nc)
        zᵢ = 0.0
        for (a, b) in zip(λ, S)
            zₖ = z[b]
            zᵢ += a * zₖ
        end
        xᵢ, yᵢ = getxy(r)
        _S = [get_neighbours(tri, i) for i in nc.indices]
        _S1 = copy(nc.indices)
        push!(_S1, reduce(union, _S)...)
        filter!(!DT.is_ghost_vertex, _S1)
        unique!(_S1)
        S = _S1
    end

    r = (xᵢ, yᵢ)
    x = [getx(get_point(tri, s)) - getx(r) for s in S]
    y = [gety(get_point(tri, s)) - gety(r) for s in S]
    X = @. [x y x^2 / 2 y^2 / 2 x * y]
    W = [1 / norm(r .- get_point(tri, s)) for s in S]
    W .= W .^ 2
    W = Diagonal(W)
    Z = [z[s] - zᵢ for s in S]
    ∇ℋ = (X' * W * X) \ (X' * W * Z)

    prob = OptimizationProblem(
        objective_function_gradient_hessian_quadratic_direct,
        copy(∇ℋ),
        (tri, S, z, xᵢ, yᵢ, zᵢ)
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD)
    return (sol.u[1:2], sol.u[3:5]), (∇ℋ[1:2], ∇ℋ[3:5])
end

function objective_function_gradient_hessian_from_initial_gradients(θ, p)
    tri, S, z, xᵢ, yᵢ, zᵢ, λ, initial_gradients, α = p
    ℓ = 0.0
    ∇₁, ∇₂, H₁₁, H₂₂, H₁₂ = θ
    for (λₛ, s) in zip(λ, S)
        zₛ = z[s]
        ∇ₛ¹, ∇ₛ² = initial_gradients[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        H = [H₁₁ H₁₂; H₁₂ H₂₂]
        ∇ = [∇₁; ∇₂]
        ∇ᵢ = [∇ₛ¹; ∇ₛ²]
        x = [xₛ - xᵢ; yₛ - yᵢ]
        L₁ = α * (zᵢ + ∇' * x + x' * H * x / 2 - zₛ)^2
        L₂ = (1 - α) * norm(H * x + ∇ - ∇ᵢ)^2
        βᵢ = λₛ / norm(x)^2
        ℓ = ℓ + βᵢ * (L₁ + L₂)
    end
    return ℓ
end
function estimate_gradient_hessian_from_initial_gradients(tri, r, z, α=0.1; use_sibson_weight=false, initial_gradients, rng=Random.default_rng())
    if r isa Integer
        S = DT.iterated_neighbourhood(tri, r, 1)
        zᵢ = z[r]
        xᵢ, yᵢ = getxy(get_point(tri, r))
        S = collect(S)
        λ = similar(S)
        fill!(λ, 1)
    else
        nc = NNI.compute_natural_coordinates(Sibson(), tri, r, rng=rng)
        S = NNI.get_indices(nc)
        λ = NNI.get_coordinates(nc)
        zᵢ = 0.0
        for (a, b) in zip(λ, S)
            zₖ = z[b]
            zᵢ += a * zₖ
        end
        xᵢ, yᵢ = getxy(r)
        if !use_sibson_weight
            fill!(λ, 1.0)
        end
    end

    r = (xᵢ, yᵢ)
    X1 = [getx(get_point(tri, s)) - getx(r) for s in S]
    X2 = [gety(get_point(tri, s)) - gety(r) for s in S]
    A = @. [X1 X2 X1^2 / 2 X2^2 / 2 X1 * X2]
    for i in axes(A, 1)
        A[i, :] .*= sqrt(α * λ[i] / norm(r .- get_point(tri, S[i]))^2)
    end
    B = [ones(size(A, 1), 1) zeros(size(A, 1), 1) X1 zeros(size(A, 1), 1) X2]
    for i in axes(B, 1)
        B[i, :] .*= sqrt((1 - α) * λ[i] / norm(r .- get_point(tri, S[i]))^2)
    end
    C = [zeros(size(A, 1), 1) ones(size(A, 1), 1) zeros(size(A, 1), 1) X2 X1]
    for i in axes(C, 1)
        C[i, :] .*= sqrt((1 - α) * λ[i] / norm(r .- get_point(tri, S[i]))^2)
    end
    D = [A; B; C]
    c = vcat(zᵢ .- z[S], [initial_gradients[s][1] for s in S], [initial_gradients[s][2] for s in S])
    for i in axes(A, 1)
        c[i] *= sqrt(α * λ[i] / norm(r .- get_point(tri, S[i]))^2)
        c[i+length(S)] *= sqrt((1 - α) * λ[i] / norm(r .- get_point(tri, S[i]))^2)
        c[i+2*length(S)] *= sqrt((1 - α) * λ[i] / norm(r .- get_point(tri, S[i]))^2)
    end
    ∇ℋ = D \ c

    prob = OptimizationProblem(
        objective_function_gradient_hessian_from_initial_gradients,
        copy(∇ℋ),
        (tri, S, z, xᵢ, yᵢ, zᵢ, λ, initial_gradients, α)
    )
    sol = solve(prob, NLopt.LN_NELDERMEAD)
    return (sol.u[1:2], sol.u[3:5]), (∇ℋ[1:2], ∇ℋ[3:5])
end

function slow_test_derivative(;
    x,
    y,
    tri,
    z,
    order,
    method,
    interpolant_method,
    alpha,
    use_cubic_terms,
    use_sibson_weight,
    initial_gradients,
    rng)
    itp = interpolate(tri, z)
    if order == 1
        ∇opt, ∇ls = estimate_gradient_direct(tri, (x, y), z; use_sibson_weight=use_sibson_weight, rng=deepcopy(rng))
        _rng = deepcopy(rng)
        S = Set{Int64}()
        S′ = Set{Int64}()
        λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, (x, y), 1, NNI.NaturalNeighboursCache(tri); rng=_rng)
        ∇man = NNI.generate_first_order_derivatives(
            method,
            tri,
            z,
            itp(x, y; method=interpolant_method, rng=_rng),
            (x, y),
            λ,
            E,
            use_sibson_weight=use_sibson_weight)
        return ∇opt, ∇ls, collect(∇man)
    else
        if method == Direct()
            if use_cubic_terms
                (∇opt, Hopt), (∇ls, Hls) = estimate_gradient_hessian_cubic_direct(tri, (x, y), z, rng=deepcopy(rng))
            else
                (∇opt, Hopt), (∇ls, Hls) = estimate_gradient_hessian_quadratic_direct(tri, (x, y), z, rng=deepcopy(rng))
            end
            _rng = deepcopy(rng)
            S = Set{Int64}()
            S′ = Set{Int64}()
            λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, (x, y), 2 + use_cubic_terms, NNI.NaturalNeighboursCache(tri); rng=_rng)
            ∇man, Hman = NNI.generate_second_order_derivatives(
                method,
                tri,
                z,
                itp(x, y; method=interpolant_method, rng=_rng),
                (x, y),
                λ,
                E,
                use_cubic_terms=use_cubic_terms
            )
            return (∇opt, Hopt), (∇ls, Hls), (collect(∇man), collect(Hman))
        else
            (∇opt, Hopt), (∇ls, Hls) = estimate_gradient_hessian_from_initial_gradients(tri, (x, y), z, alpha; use_sibson_weight=use_sibson_weight, initial_gradients=initial_gradients, rng=deepcopy(rng))
            _rng = deepcopy(rng)
            S = Set{Int64}()
            S′ = Set{Int64}()
            λ, E = NNI.get_taylor_neighbourhood!(S, S′, tri, (x, y), 1, NNI.NaturalNeighboursCache(tri); rng=_rng)
            ∇man, Hman = NNI.generate_second_order_derivatives(
                method,
                tri,
                z,
                itp(x, y; method=interpolant_method, rng=_rng),
                (x, y),
                λ,
                E,
                alpha=alpha,
                use_cubic_terms=use_cubic_terms,
                initial_gradients=initial_gradients
            )
            return (∇opt, Hopt), (∇ls, Hls), (collect(∇man), collect(Hman))
        end
    end
    throw("Invalid arguments.")
end

function slow_test_derivative(itp1, itp2, ∂11, ∂12, ∂21, ∂22; x, y, rng, tri, z, kwargs...)
    grad11 = collect(∂11(x, y; rng=deepcopy(rng), kwargs...))
    grad12, hess12 = collect.(∂12(x, y; rng=deepcopy(rng), kwargs...))
    grad21 = collect(∂21(x, y; rng=deepcopy(rng), kwargs...))
    grad22, hess22 = collect.(∂22(x, y; rng=deepcopy(rng), kwargs...))

    gradopt11, gradls11, gradman11 = collect.(slow_test_derivative(; x, y, rng=deepcopy(rng), tri, z, order=1, initial_gradients=itp1.gradient, kwargs...))
    (gradopt12, hessopt12), (gradls12, hessls12), (gradman12, hessman12) = collect.(slow_test_derivative(; x, y, tri, z, rng=deepcopy(rng), order=2, initial_gradients=itp1.gradient, kwargs...))
    gradopt21, gradls21, gradman21 = collect.(slow_test_derivative(; x, y, rng=deepcopy(rng), order=1, tri, z, initial_gradients=itp2.gradient, kwargs...))
    (gradopt22, hessopt22), (gradls22, hessls22), (gradman22, hessman22) = collect.(slow_test_derivative(; x, y, tri, z, rng=deepcopy(rng), order=2, initial_gradients=itp2.gradient, kwargs...))

    flag1 = isapprox(grad11, gradman11, rtol=1e-1)
    flag2 = isapprox(grad12, gradman12, rtol=1e-1)
    flag3 = isapprox(hess12, hessman12, rtol=1e-1)
    flag2 = isapprox(grad12, gradman12, rtol=1e-1)
    flag3 = isapprox(hess12, hessman12, rtol=1e-1)
    flag4 = isapprox(grad21, gradman21, rtol=1e-1)
    flag5 = isapprox(grad22, gradman22, rtol=1e-1)
    flag6 = isapprox(hess22, hessman22, rtol=1e-1)
    flag7 = isapprox(grad11, gradopt11, rtol=1e-1)
    flag8 = isapprox(grad12, gradopt12, rtol=1e-1)
    flag9 = isapprox(hess12, hessopt12, rtol=1e-1)
    flag10 = isapprox(grad21, gradopt21, rtol=1e-1)
    flag11 = isapprox(grad22, gradopt22, rtol=1e-1)
    flag12 = isapprox(hess22, hessopt22, rtol=1e-1)
    flag13 = isapprox(grad11, gradls11, rtol=1e-1)
    flag14 = isapprox(grad12, gradls12, rtol=1e-1)
    flag15 = isapprox(hess12, hessls12, rtol=1e-1)
    flag16 = isapprox(grad21, gradls21, rtol=1e-1)
    flag17 = isapprox(grad22, gradls22, rtol=1e-1)
    flag18 = isapprox(hess22, hessls22, rtol=1e-1)
    flag19 = ∂11(x, y; rng=deepcopy(rng), kwargs...) |> collect == ∂11(x, y; rng=deepcopy(rng), kwargs...) |> collect
    flag20 = ∂12(x, y; rng=deepcopy(rng), kwargs...) |> collect == ∂12(x, y; rng=deepcopy(rng), kwargs...) |> collect
    flag21 = ∂21(x, y; rng=deepcopy(rng), kwargs...) |> collect == ∂21(x, y; rng=deepcopy(rng), kwargs...) |> collect
    flag22 = ∂22(x, y; rng=deepcopy(rng), kwargs...) |> collect == ∂22(x, y; rng=deepcopy(rng), kwargs...) |> collect
    all_flags = (flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8, flag9,
        flag10, flag11, flag12, flag13, flag14, flag15, flag16, flag17, flag18,
        flag19, flag20, flag21, flag22)
    final_flag = all(all_flags)
    if !final_flag
        idx = findall(!, all_flags)
        println("Failed at index: $idx")
        return false
    end
    return true
end