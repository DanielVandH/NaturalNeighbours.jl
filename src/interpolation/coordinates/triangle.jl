function _compute_triangle_coordinates(
    tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL},
    interpolation_point,
    cache::NaturalNeighboursCache{F}=NaturalNeighboursCache(tri);
    project=true,
    kwargs...
) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,F}
    coordinates = get_coordinates(cache)
    envelope = get_envelope(cache)
    last_triangle = get_last_triangle(cache)
    V = jump_and_march(tri, interpolation_point; try_points=last_triangle[], kwargs...)
    i, j, return_flag = check_for_extrapolation(tri, V, interpolation_point, last_triangle)
    return_flag && return two_point_interpolate!(coordinates, envelope, tri, i, j, interpolation_point, project)
    i, j, k = indices(V)
    resize!(coordinates, 3)
    resize!(envelope, 3)
    p, q, r = get_point(tri, i, j, k)
    x₁, y₁ = getxy(p)
    x₂, y₂ = getxy(q)
    x₃, y₃ = getxy(r)
    x, y = getxy(interpolation_point)
    Δ = (y₂ - y₃) * (x₁ - x₃) + (x₃ - x₂) * (y₁ - y₃)
    λ₁ = ((y₂ - y₃) * (x - x₃) + (x₃ - x₂) * (y - y₃)) / Δ
    λ₂ = ((y₃ - y₁) * (x - x₃) + (x₁ - x₃) * (y - y₃)) / Δ
    λ₃ = one(λ₁) - λ₁ - λ₂
    coordinates[1] = λ₁
    coordinates[2] = λ₂
    coordinates[3] = λ₃
    envelope[1] = i
    envelope[2] = j
    envelope[3] = k
    return NaturalCoordinates(coordinates, envelope, interpolation_point, tri)
end

function compute_natural_coordinates(::Triangle, tri, interpolation_point, cache=NaturalNeighboursCache(tri); kwargs...) 
    return _compute_triangle_coordinates(tri, interpolation_point, cache; kwargs...)
end