function two_point_interpolate!(tri, i, j, c)
    # Project c onto the line through (a, b): https://stackoverflow.com/a/15187473
    # The orthogonal projection is not guaranteed to be on the line segment (corresponding 
    # to t < 0 or t > 1), in which case the weights are no longer a convex combination.
    a, b = get_point(tri, i, j)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    ℓ² = (ax - bx)^2 + (ay - by)^2
    t = (cx - ax) * (bx - ax) + (cy - ay) * (by - ay)
    t /= ℓ²
    return t
end

function two_point_interpolate!(coordinates::AbstractVector{F}, envelope, tri, i, j, r, project=true) where {F} #interpolate r using two points i, j
    if project
        t = two_point_interpolate!(tri, i, j, r)
        resize!(coordinates, 2)
        resize!(envelope, 2)
        coordinates[1] = one(t) - t
        coordinates[2] = t
        envelope[1] = i
        envelope[2] = j
    else
        resize!(coordinates, 1)
        resize!(envelope, 1)
        coordinates[1] = F(NaN)
        envelope[1] = i
    end
    return NaturalCoordinates(coordinates, envelope, r, tri)
end

function check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    if is_ghost_triangle(V)
        V = rotate_ghost_triangle_to_standard_form(V)
        i, j, _ = indices(V)
        last_triangle[] = (i, j, get_adjacent(tri, j, i))
    else
        last_triangle[] = indices(V)
    end
    if is_boundary_triangle(tri, V)
        _V = replace_boundary_triangle_with_ghost_triangle(tri, V)
        _u, _w, _ = indices(_V)
        cert = point_position_relative_to_line(tri, _u, _w, interpolation_point)
        if is_collinear(cert)
            cert = point_position_on_line_segment(tri, _u, _w, interpolation_point)
            if is_on(cert) || is_degenerate(cert)
                V = _V
            end
        end
    end
    if is_ghost_triangle(V)
        V = rotate_ghost_triangle_to_standard_form(V)
        i, j, _ = indices(V)
        return i, j, true
    end
    i, j, _ = indices(V)
    return i, j, false
end
