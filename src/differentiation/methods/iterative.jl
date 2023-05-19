function _generate_second_order_derivatives_iterative(
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    initial_gradients,
    d_cache=DerivativeCache(tri),
    alpha=0.1,
    use_sibson_weight=true
)
    X = get_quadratic_matrix_no_cubic(d_cache)
    b = get_rhs_vector(d_cache)
    ∇ℋ = get_quadratic_sol_no_cubic(d_cache)
    p = get_point(tri, i)
    xᵢ, yᵢ = getxy(p)
    m = length(E)
    resize!(X, 5, 3m)
    resize!(b, 3m)
    α = alpha
    α′ = one(alpha) - alpha
    for (j, s) in enumerate(E)
        λₛ = get_λ(λ, j, use_sibson_weight)
        zₛ = z[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        ∇ₛ¹, ∇ₛ² = initial_gradients[s]
        δ = (xᵢ - xₛ)^2 + (yᵢ - yₛ)^2
        βₛ = λₛ * inv(δ)
        γₛ = sqrt(α * βₛ)
        γₛ′ = sqrt(α′ * βₛ)
        X[1, j] = γₛ * (xₛ - xᵢ)
        X[2, j] = γₛ * (yₛ - yᵢ)
        X[3, j] = γₛ * (xₛ - xᵢ)^2 / 2
        X[4, j] = γₛ * (yₛ - yᵢ)^2 / 2
        X[5, j] = γₛ * (xₛ - xᵢ) * (yₛ - yᵢ)
        j′ = m + j
        X[1, j′] = γₛ′
        X[2, j′] = zero(γₛ′)
        X[3, j′] = γₛ′ * (xₛ - xᵢ)
        X[4, j′] = zero(γₛ′)
        X[5, j′] = γₛ′ * (yₛ - yᵢ)
        j′′ = 2m + j
        X[1, j′′] = zero(γₛ′)
        X[2, j′′] = γₛ′
        X[3, j′′] = zero(γₛ′)
        X[4, j′′] = γₛ′ * (yₛ - yᵢ)
        X[5, j′′] = γₛ′ * (xₛ - xᵢ)
        z̃ₛ = zₛ - zᵢ
        b[j] = γₛ * z̃ₛ
        b[j′] = γₛ′ * ∇ₛ¹
        b[j′′] = γₛ′ * ∇ₛ²
    end
    qr_X = qr!(X')
    ldiv!(∇ℋ, qr_X, b)
    return (∇ℋ[1], ∇ℋ[2]), (∇ℋ[3], ∇ℋ[4], ∇ℋ[5])
end

function generate_second_order_derivatives(
    method::Iterative,
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri);
    alpha=0.1,
    initial_gradients=nothing,
    use_cubic_terms=nothing, # not used,
    use_sibson_weight=true
)
    if isnothing(initial_gradients)
        throw(ArgumentError("initial_gradients must be provided for iterative derivative estimation. Consider using e.g. interpolate(tri, z; derivatives = true)."))
    end
    return _generate_second_order_derivatives_iterative(tri, z, zᵢ, i, λ, E, initial_gradients, d_cache, alpha, use_sibson_weight)
end