function _eval_gradient_direct(
    tri,
    z,
    i,
    cache=DerivativeCache(tri);
)
    S = get_iterated_neighbourhood(cache)
    A = get_linear_matrix(cache)
    b = get_rhs_vector(cache)
    ∇ = get_linear_sol(cache)
    p = get_point(tri, i)
    x, y = getxy(p)
    iterated_neighbourhood!(S, tri, i, 1)
    m = length(S)
    resize!(A, (3, m))
    resize!(b, m)
    for (i, s) in enumerate(S)
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        δ² = (x - xₛ)^2 + (y - yₛ)^2
        δ = sqrt(δ²)
        βᵢ = inv(δ)
        aᵢ₁ = βᵢ
        aᵢ₂ = βᵢ * xₛ
        aᵢ₃ = βᵢ * yₛ
        A[1, i] = aᵢ₁
        A[2, i] = aᵢ₂
        A[3, i] = aᵢ₃
        b[i] = βᵢ * z[s]
    end
    qr_A = qr!(A')
    ldiv!(∇, qr_A, b)
    return (∇[2], ∇[3])
end

# gradients are hard to estimate for some reason when 
# we also use a Hessian, but they are ok when we include 
# cubic terms - consistent with the discussion 
# on p. 100 of Bobach's thesis. So, we do waste some 
# information here computing cubic terms, but it's
# not necessary and not a big deal anyway.
function _eval_gradient_and_hessian_direct(
    tri,
    z,
    i,
    cache=DerivativeCache(tri);
    kwargs...
)
    S = get_iterated_neighbourhood(cache)
    A = get_quadratic_matrix(cache)
    b = get_rhs_vector(cache)
    ∇ℋ = get_quadratic_sol(cache)
    p = get_point(tri, i)
    x, y = getxy(p)
    iterated_neighbourhood!(S, tri, i, 3)
    m = length(S)
    resize!(A, (10, m))
    resize!(b, m)
    for (j, s) in enumerate(S)
        pₛ = get_point(tri, s)
        xₛ, yₛ = getxy(pₛ)
        δ² = (x - xₛ)^2 + (y - yₛ)^2
        δ = sqrt(δ²)
        βᵢ = inv(δ)
        aᵢ₁ = βᵢ
        aᵢ₂ = βᵢ * xₛ
        aᵢ₃ = βᵢ * yₛ
        aᵢ₄ = βᵢ * xₛ^2 / 2
        aᵢ₅ = βᵢ * yₛ^2 / 2
        aᵢ₆ = βᵢ * xₛ * yₛ
        aᵢ₇ = βᵢ * xₛ^3 / 6
        aᵢ₈ = βᵢ * yₛ^3 / 6
        aᵢ₉ = βᵢ * xₛ^2 * yₛ / 2
        aᵢ₁₀ = βᵢ * xₛ * yₛ^2 / 2
        A[1, j] = aᵢ₁
        A[2, j] = aᵢ₂
        A[3, j] = aᵢ₃
        A[4, j] = aᵢ₄
        A[5, j] = aᵢ₅
        A[6, j] = aᵢ₆
        A[7, j] = aᵢ₇
        A[8, j] = aᵢ₈
        A[9, j] = aᵢ₉
        A[10, j] = aᵢ₁₀
        b[j] = βᵢ * (z[s] - z[i])
    end
    qr_A = qr!(A')
    ldiv!(∇ℋ, qr_A, b)
    return (∇ℋ[2], ∇ℋ[3]), (∇ℋ[4], ∇ℋ[5], ∇ℋ[6])
end

function eval_gradient(::Direct, tri::Triangulation, z, p, cache=DerivativeCache(tri))
    return _eval_gradient_direct(tri, z, p, cache)
end

function eval_gradient_and_hessian(::Direct, tri::Triangulation, z, p, cache=DerivativeCache(tri))
    return _eval_gradient_and_hessian_direct(tri, z, p, cache)
end