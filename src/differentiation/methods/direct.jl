function _generate_first_order_derivatives_direct(
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri, z);
    use_sibson_weight=!(i isa Integer))
    X = get_linear_matrix(d_cache)
    b = get_rhs_vector(d_cache)
    p = get_point(tri, i)
    xᵢ, yᵢ = getxy(p)
    m = length(E)
    resize!(X, 2, m)
    is_scalar(z) ? resize!(b, m) : resize!(b, fdim(z), m)
    for (j, s) in enumerate(E)
        λₛ = get_λ(λ, j, use_sibson_weight)
        zₛ = get_data(z, s)
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        δ = sqrt((xᵢ - xₛ)^2 + (yᵢ - yₛ)^2)
        βₛ = sqrt(λₛ) * inv(δ)
        X[1, j] = βₛ * (xₛ - xᵢ)
        X[2, j] = βₛ * (yₛ - yᵢ)
        if is_scalar(z)
            b[j] = βₛ * (zₛ - zᵢ)
        else
            @. b[:, j] = βₛ * (zₛ - zᵢ)
        end
    end
    #=
    @static if VERSION < v"1.7.0"
        qr_X = qr!(Matrix(X'))
    else
        qr_X = qr!(X')
    end
    =#
    ∇ = is_scalar(z) ? X' \ b : X' \ b' # I want to do the above, but it just keeps segfaulting somehow
    if is_scalar(z)
        return (∇[1], ∇[2])
    else
        gradient = get_gradient(d_cache)
        @. gradient = ∇[1:2, :]
        return gradient
    end
end

function _generate_second_order_derivatives_direct(
    tri,
    z,
    zᵢ,
    i,
    E,
    d_cache=DerivativeCache(tri, z);
    use_cubic_terms=true)
    X = if use_cubic_terms
        get_quadratic_matrix(d_cache)
    else
        get_quadratic_matrix_no_cubic(d_cache)
    end
    b = get_rhs_vector(d_cache)
    p = get_point(tri, i)
    xᵢ, yᵢ = getxy(p)
    m = length(E)
    resize!(X, 5 + 4use_cubic_terms, m)
    is_scalar(z) ? resize!(b, m) : resize!(b, fdim(z), m)
    for (j, s) in enumerate(E)
        zₛ = get_data(z, s)
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        δ = sqrt((xᵢ - xₛ)^2 + (yᵢ - yₛ)^2)
        βₛ = inv(δ)
        X[1, j] = βₛ * (xₛ - xᵢ)
        X[2, j] = βₛ * (yₛ - yᵢ)
        X[3, j] = βₛ * (xₛ - xᵢ)^2 / 2
        X[4, j] = βₛ * (yₛ - yᵢ)^2 / 2
        X[5, j] = βₛ * (xₛ - xᵢ) * (yₛ - yᵢ)
        if use_cubic_terms
            X[6, j] = βₛ * (xₛ - xᵢ)^3 / 6
            X[7, j] = βₛ * (yₛ - yᵢ)^3 / 6
            X[8, j] = βₛ * (xₛ - xᵢ)^2 * (yₛ - yᵢ) / 2
            X[9, j] = βₛ * (xₛ - xᵢ) * (yₛ - yᵢ)^2 / 2
        end
        if is_scalar(z)
            b[j] = βₛ * (zₛ - zᵢ)
        else
            @. b[:, j] = βₛ * (zₛ - zᵢ)
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

function generate_first_order_derivatives(
    method::AbstractDifferentiator, # Iterative() also goes here
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri, z);
    use_sibson_weight=!(i isa Integer),
    kwargs...
)
    return _generate_first_order_derivatives_direct(tri, z, zᵢ, i, λ, E, d_cache; use_sibson_weight)
end

function generate_second_order_derivatives(
    ::Direct,
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri, z);
    use_cubic_terms=true,
    kwargs...
)
    return _generate_second_order_derivatives_direct(tri, z, zᵢ, i, E, d_cache; use_cubic_terms)
end


