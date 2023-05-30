function _generate_first_order_derivatives_direct(
    tri,
    z,
    zᵢ,
    i,
    λ,
    E,
    d_cache=DerivativeCache(tri);
    use_sibson_weight=!(i isa Integer))
    X = get_linear_matrix(d_cache)
    b = get_rhs_vector(d_cache)
    p = get_point(tri, i)
    xᵢ, yᵢ = getxy(p)
    m = length(E)
    resize!(X, 2, m)
    resize!(b, m)
    for (j, s) in enumerate(E)
        λₛ = get_λ(λ, j, use_sibson_weight)
        zₛ = z[s]
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        δ = sqrt((xᵢ - xₛ)^2 + (yᵢ - yₛ)^2)
        βₛ = sqrt(λₛ) * inv(δ)
        X[1, j] = βₛ * (xₛ - xᵢ)
        X[2, j] = βₛ * (yₛ - yᵢ)
        b[j] = βₛ * (zₛ - zᵢ)
    end
    @static if VERSION < v"1.7.0"
        qr_X = qr!(Matrix(X'))
    else
        qr_X = qr!(X')
    end
    ∇ = copy(b) # This is the same fix in https://github.com/JuliaLang/julia/pull/43510 to avoid views, avoiding shared data issues
    ldiv!(qr_X, ∇)
    return (∇[1], ∇[2])
end

function _generate_second_order_derivatives_direct(
    tri,
    z,
    zᵢ,
    i,
    E,
    d_cache=DerivativeCache(tri);
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
    resize!(b, m)
    for (j, s) in enumerate(E)
        zₛ = z[s]
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
        b[j] = βₛ * (zₛ - zᵢ)
    end
    @static if VERSION < v"1.7.0"
        ∇ℋ = X' \ b
        return (∇ℋ[1], ∇ℋ[2]), (∇ℋ[3], ∇ℋ[4], ∇ℋ[5])
    else
        qr_X = qr!(X')
        ∇ℋ = copy(b) # This is the same fix in https://github.com/JuliaLang/julia/pull/43510 to avoid views, avoiding shared data issues
        5 + 4use_cubic_terms > m && resize!(∇ℋ, 5 + 4use_cubic_terms) # See Issue #11
        ldiv!(qr_X, ∇ℋ)
        return (∇ℋ[1], ∇ℋ[2]), (∇ℋ[3], ∇ℋ[4], ∇ℋ[5])
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
    d_cache=DerivativeCache(tri);
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
    d_cache=DerivativeCache(tri);
    use_cubic_terms=true,
    kwargs...
)
    return _generate_second_order_derivatives_direct(tri, z, zᵢ, i, E, d_cache; use_cubic_terms)
end


