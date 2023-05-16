## gradient 
function prepare_gradient_evaluation(tri, i, cache=DerivativeCache(tri); use_sibson_weight=false)
    S = get_iterated_neighbourhood(cache)
    X = get_linear_matrix(cache)
    b = get_rhs_vector(cache)
    ∇ = get_linear_sol(cache)
    p = get_point(tri, i)
    x, y = getxy(p)
    iterated_neighbourhood!(S, tri, i, 1)
    m = length(S)
    resize!(X, (2, m))
    resize!(b, m)
    return S, X, b, ∇, x, y
end

function get_residual_weight(tri, s, x, y, λ=nothing)
    pₛ = get_point(tri, s)
    xₛ, yₛ = getxy(pₛ)
    δ² = (x - xₛ)^2 + (y - yₛ)^2
    βₛ = inv(δ²)
    if !isnothing(λ)
        βₛ = βₛ * λ
    end
    return βₛ, xₛ, yₛ
end

function fill_gradient_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
    X[1, j] = βₛ * (xₛ - x)
    X[2, j] = βₛ * (yₛ - y)
    return nothing
end

function solve_gradient_problem!(X, b, ∇)
    qr_X = qr!(X')
    ldiv!(∇, qr_X, b)
    return (∇[1], ∇[2])
end

function _eval_gradient_direct(
    tri,
    z,
    i,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false # not currently used. No idea what λ₁ˢⁱᵇ⁰ is in the thesis - it's just 1?! 
)
    S, X, b, ∇, x, y = prepare_gradient_evaluation(tri, i, d_cache)
    for (j, s) in enumerate(S)
        βₛ, xₛ, yₛ = get_residual_weight(tri, s, x, y)
        βₛ = sqrt(βₛ)
        fill_gradient_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
        b[j] = βₛ * (z[s] - z[i])
    end
    return solve_gradient_problem!(X, b, ∇)
end

function eval_gradient(
    ::Direct,
    tri::Triangulation,
    z,
    p,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false, kwargs...)
    return _eval_gradient_direct(tri, z, p, d_cache, n_cache; use_sibson_weight, kwargs...)
end

## hessian
function prepare_hessian_evaluation(tri, i, cache=DerivativeCache(tri); use_cubic_terms=true)
    S = get_iterated_neighbourhood(cache)
    X = get_quadratic_matrix(cache)
    b = get_rhs_vector(cache)
    ∇ℋ = get_quadratic_sol(cache)
    p = get_point(tri, i)
    x, y = getxy(p)
    iterated_neighbourhood!(S, tri, i, 2 + use_cubic_terms) # 2 or 3
    m = length(S)
    resize!(X, (9, m))
    resize!(b, m)
    return S, X, b, ∇ℋ, x, y
end

function fill_hessian_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
    X[3, j] = βₛ * (xₛ - x)^2 / 2
    X[4, j] = βₛ * (yₛ - y)^2 / 2
    X[5, j] = βₛ * (xₛ - x) * (yₛ - y)
    return nothing
end

function fill_cubic_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
    X[6, j] = βₛ * (xₛ - x)^3 / 6
    X[7, j] = βₛ * (yₛ - y)^3 / 6
    X[8, j] = βₛ * (xₛ - x)^2 * (yₛ - y) / 2
    X[9, j] = βₛ * (xₛ - x) * (yₛ - y)^2 / 2
    return nothing
end

function solve_gradient_hessian_problem!(X, b, ∇ℋ; use_cubic_terms=true)
    if !use_cubic_terms
        @views X = X[1:5, :]
    end
    X = @views X[1:(5+4use_cubic_terms), :] # 5 or 9
    qr_X = qr!(X')
    ldiv!(∇ℋ, qr_X, b)
    return (∇ℋ[1], ∇ℋ[2]), (∇ℋ[3], ∇ℋ[4], ∇ℋ[5])
end

# gradients are hard to estimate for some reason when 
# we also use a Hessian, but they are ok when we include 
# cubic terms - consistent with the discussion 
# on p. 100 of Bobach's thesis. That's why we 
# use use_cubic_terms by default. Disable it if you please.
function _eval_gradient_and_hessian_direct(
    tri,
    z,
    i,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false,
    use_cubic_terms=true,
    kwargs...
)
    S, X, b, ∇ℋ, x, y = prepare_hessian_evaluation(tri, i, d_cache; use_cubic_terms)
    for (j, s) in enumerate(S)
        βₛ, xₛ, yₛ = get_residual_weight(tri, s, x, y)
        βₛ = sqrt(βₛ)
        fill_gradient_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
        fill_hessian_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
        use_cubic_terms && fill_cubic_coefficients!(X, βₛ, x, y, xₛ, yₛ, j)
        b[j] = βₛ * (z[s] - z[i])
    end
    return solve_gradient_hessian_problem!(X, b, ∇ℋ; use_cubic_terms)
end

function eval_gradient_and_hessian(
    ::Direct,
    tri::Triangulation,
    z,
    p,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false,
    use_cubic_terms=true,
    kwargs...)
    return _eval_gradient_and_hessian_direct(tri, z, p, d_cache, n_cache; use_sibson_weight, use_cubic_terms, kwargs...)
end