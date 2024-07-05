function _compute_triangle_shape_coefficients(tri, i, j, k)
    p, q, r = get_point(tri, i, j, k)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    Δ = qx * ry - qy * rx - px * ry + rx * py + px * qy - qx * py
    s₁ = (qy - ry) / Δ
    s₂ = (ry - py) / Δ
    s₃ = (py - qy) / Δ
    s₄ = (rx - qx) / Δ
    s₅ = (px - rx) / Δ
    s₆ = (qx - px) / Δ
    s₇ = (qx * ry - rx * qy) / Δ
    s₈ = (rx * py - px * ry) / Δ
    s₉ = (px * qy - qx * py) / Δ
    shape_function_coefficients = (s₁, s₂, s₃, s₄, s₅, s₆, s₇, s₈, s₉)
    return shape_function_coefficients
end