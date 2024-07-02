function _generate_second_order_derivatives_iterative(
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    initial_gradients,
    d_cache=DerivativeCache(tri, z),
    alpha=0.1,
    use_sibson_weight=true
)
    X = get_quadratic_matrix_no_cubic(d_cache)
    b = get_rhs_vector(d_cache)
    p = get_point(tri, i)
    xᵢ, yᵢ = getxy(p)
    m = length(E)
    resize!(X, 5, 3m)
    is_scalar(z) ? resize!(b, 3m) : resize!(bm, fdim(z), 3m)
    α = alpha
    α′ = one(alpha) - alpha
    for (j, s) in enumerate(E)
        λₛ = get_λ(λ, j, use_sibson_weight)
        zₛ = get_data(z, s)
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        ∇ₛ¹² = get_data(initial_gradients, s)
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
        if is_scalar(z)
            z̃ₛ = zₛ - zᵢ
            b[j] = γₛ * z̃ₛ
            b[j′] = γₛ′ * ∇ₛ¹²[1]
            b[j′′] = γₛ′ * ∇ₛ¹²[2]
        else
            @. b[:, j] = γₛ * (zₛ - zᵢ)
            for idx in 1:fdim(z)
                b[j′, idx] = γₛ′ * ∇ₛ¹²[idx][1]
                b[j′′, idx] = γₛ′ * ∇ₛ¹²[idx][2]
            end
        end
    end
    if is_scalar(z)
        ∇ℋ = X' \ b
        return (∇ℋ[1], ∇ℋ[2]), (∇ℋ[3], ∇ℋ[4], ∇ℋ[5])
    else
        ∇ℋ = X' \ b'
        gradient = get_gradient(d_cache)
        hessian = get_hessian(d_cache)
        @. gradient = ∇ℋ[1:2, :]
        @. hessian = ∇ℋ[3:5, :]
        return gradient, hessian
    end
end

function generate_second_order_derivatives(
    method::Iterative,
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri, z);
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