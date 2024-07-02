function _compute_hiyoshi_coordinates(nc::NaturalCoordinates{F}, tri::Triangulation, z, ∇, H) where {F}
    λ = get_coordinates(nc)
    N₀ = get_indices(nc)
    result = is_scalar(z) ? zero(F) : zeros(F, fdim(z))
    for i in eachindex(λ)
        for j in i:lastindex(λ)
            for k in j:lastindex(λ)
                for ℓ in k:lastindex(λ)
                    for m in ℓ:lastindex(λ)
                        (i′, j′, k′, ℓ′, m′), case = group_sort(i, j, k, ℓ, m) # could (??) be faster to sort individually at each loop step, but not bothered with that currently
                        bezier_height, multiplicity = get_contrib(tri, i′, j′, k′, ℓ′, m′, N₀, ∇, H, z, case)
                        λ_prod = λ[i′] * λ[j′] * λ[k′] * λ[ℓ′] * λ[m′]
                        if is_scalar(z)
                            result += 120bezier_height * λ_prod / multiplicity
                        else
                            @. result += 120bezier_height * λ_prod / multiplicity
                        end
                    end
                end
            end
        end
    end
    return result
end

function _hiyoshi_case_1(i, N₀, z) # iiiii
    zᵢ = get_data(z, N₀[i])
    fᵢᵢᵢᵢᵢ = zᵢ
    return fᵢᵢᵢᵢᵢ
end
function _hiyoshi_case_2(tri, i, j, N₀, ∇, z) # iiiij
    zᵢ = get_data(z, N₀[i])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    fᵢᵢᵢᵢⱼ = @. zᵢ + zᵢⱼ / 5
    return fᵢᵢᵢᵢⱼ
end
function _hiyoshi_case_3(tri, i, j, N₀, ∇, H, z) # iiijj
    zᵢ = get_data(z, N₀[i])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    zᵢⱼⱼ = hessian_form(tri, i, j, j, N₀, H)
    fᵢᵢᵢⱼⱼ = @. zᵢ + 2zᵢⱼ / 5 + zᵢⱼⱼ / 20
    return fᵢᵢᵢⱼⱼ
end
function _hiyoshi_case_4(tri, i, j, k, N₀, ∇, H, z) # iiijk
    zᵢ = get_data(z, N₀[i])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    zᵢₖ = directional_derivative(tri, i, k, N₀, ∇)
    zᵢⱼₖ = hessian_form(tri, i, j, k, N₀, H)
    fᵢᵢᵢⱼₖ = @. zᵢ + (zᵢⱼ + zᵢₖ) / 5 + zᵢⱼₖ / 20
    return fᵢᵢᵢⱼₖ
end
function _hiyoshi_case_5(tri, i, j, k, N₀, ∇, H, z) # iijjk
    zᵢ = get_data(z, N₀[i])
    zⱼ = get_data(z, N₀[j])
    zₖ = get_data(z, N₀[k])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    zⱼᵢ = directional_derivative(tri, j, i, N₀, ∇)
    zᵢₖ = directional_derivative(tri, i, k, N₀, ∇)
    zⱼₖ = directional_derivative(tri, j, k, N₀, ∇)
    zₖᵢ = directional_derivative(tri, k, i, N₀, ∇)
    zₖⱼ = directional_derivative(tri, k, j, N₀, ∇)
    zᵢⱼₖ = hessian_form(tri, i, j, k, N₀, H)
    zⱼᵢₖ = hessian_form(tri, j, i, k, N₀, H)
    zₖᵢⱼ = hessian_form(tri, k, i, j, N₀, H)
    fᵢᵢⱼⱼₖ = @. 13(zᵢ + zⱼ) / 30 + 2zₖ / 15 + (zᵢⱼ + zⱼᵢ) / 9 + 7(zᵢₖ + zⱼₖ) / 90 + 2(zₖᵢ + zₖⱼ) / 45 + (zᵢⱼₖ + zⱼᵢₖ + zₖᵢⱼ) / 45
    # fᵢᵢⱼⱼₖ = (zᵢ + zⱼ) / 2 + 3(zᵢⱼ + zⱼᵢ) / 20 + (zᵢₖ + zⱼₖ) / 10 + (zᵢⱼₖ + zⱼᵢₖ) / 30 + (zᵢⱼⱼ + zⱼᵢᵢ) / 120
    return fᵢᵢⱼⱼₖ
end
function _hiyoshi_case_6(tri, i, j, k, ℓ, N₀, ∇, H, z) # iijkℓ
    zᵢ = get_data(z, N₀[i])
    zⱼ = get_data(z, N₀[j])
    zₖ = get_data(z, N₀[k])
    zₗ = get_data(z, N₀[ℓ])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    zᵢₖ = directional_derivative(tri, i, k, N₀, ∇)
    zᵢₗ = directional_derivative(tri, i, ℓ, N₀, ∇)
    zᵢⱼₖ = hessian_form(tri, i, j, k, N₀, H)
    zᵢⱼₗ = hessian_form(tri, i, j, ℓ, N₀, H)
    zᵢₖₗ = hessian_form(tri, i, k, ℓ, N₀, H)
    zⱼᵢ = directional_derivative(tri, j, i, N₀, ∇)
    zⱼₖ = directional_derivative(tri, j, k, N₀, ∇)
    zⱼₗ = directional_derivative(tri, j, ℓ, N₀, ∇)
    zₖᵢ = directional_derivative(tri, k, i, N₀, ∇)
    zₖⱼ = directional_derivative(tri, k, j, N₀, ∇)
    zₖₗ = directional_derivative(tri, k, ℓ, N₀, ∇)
    zₗᵢ = directional_derivative(tri, ℓ, i, N₀, ∇)
    zₗⱼ = directional_derivative(tri, ℓ, j, N₀, ∇)
    zₗₖ = directional_derivative(tri, ℓ, k, N₀, ∇)
    zⱼᵢₖ = hessian_form(tri, j, i, k, N₀, H)
    zⱼᵢₗ = hessian_form(tri, j, i, ℓ, N₀, H)
    zⱼₖₗ = hessian_form(tri, j, k, ℓ, N₀, H)
    zₖᵢⱼ = hessian_form(tri, k, i, j, N₀, H)
    zₖᵢₗ = hessian_form(tri, k, i, ℓ, N₀, H)
    zₖⱼₗ = hessian_form(tri, k, j, ℓ, N₀, H)
    zₗᵢⱼ = hessian_form(tri, ℓ, i, j, N₀, H)
    zₗᵢₖ = hessian_form(tri, ℓ, i, k, N₀, H)
    zₗⱼₖ = hessian_form(tri, ℓ, j, k, N₀, H)
    #=
    fᵢᵢⱼₖₗ = 7(zᵢ + zⱼ + zₖ + zₗ) / 10 +
             11(zᵢⱼ + zᵢₖ + zᵢₗ) / 90 +
             (zᵢⱼₖ + zᵢⱼₗ + zᵢₖₗ) / 45 +
             (zⱼᵢ + zⱼₖ + zⱼₗ + zₖᵢ + zₖⱼ + zₖₗ + zₗᵢ + zₗⱼ + zₗₖ) / 45 +
             (zⱼᵢₖ + zⱼᵢₗ + zⱼₖₗ + zₖᵢⱼ + zₖᵢₗ + zₖⱼₗ + zₗᵢⱼ + zₗᵢₖ + zₗⱼₖ) / 180
    =#
    fᵢᵢⱼₖₗ = @. zᵢ / 2 + (zⱼ + zₖ + zₗ) / 6 + 7(zᵢⱼ + zᵢₖ + zᵢₗ) / 90 + 2(zⱼᵢ + zₖᵢ + zₗᵢ) / 45 +
                (zⱼₖ + zⱼₗ + zₖⱼ + zₖₗ + zₗⱼ + zₗₖ) / 30 + (zᵢⱼₖ + zᵢⱼₗ + zᵢₖₗ) / 90 +
                (zⱼᵢₖ + zⱼᵢₗ + zₖᵢⱼ + zₖᵢₗ + zₗᵢⱼ + zₗᵢₖ) / 90 + (zⱼₖₗ + zₖⱼₗ + zₗⱼₖ) / 180
    return fᵢᵢⱼₖₗ
end
function _hiyoshi_case_7(tri, i, j, k, ℓ, m, N₀, ∇, H, z) # ijkℓm
    zᵢ = get_data(z, N₀[i])
    zⱼ = get_data(z, N₀[j])
    zₖ = get_data(z, N₀[k])
    zₗ = get_data(z, N₀[ℓ])
    zₘ = get_data(z, N₀[m])
    zᵢⱼ = directional_derivative(tri, i, j, N₀, ∇)
    zᵢₖ = directional_derivative(tri, i, k, N₀, ∇)
    zᵢₗ = directional_derivative(tri, i, ℓ, N₀, ∇)
    zᵢₘ = directional_derivative(tri, i, m, N₀, ∇)
    zⱼᵢ = directional_derivative(tri, j, i, N₀, ∇)
    zⱼₖ = directional_derivative(tri, j, k, N₀, ∇)
    zⱼₗ = directional_derivative(tri, j, ℓ, N₀, ∇)
    zⱼₘ = directional_derivative(tri, j, m, N₀, ∇)
    zₖᵢ = directional_derivative(tri, k, i, N₀, ∇)
    zₖⱼ = directional_derivative(tri, k, j, N₀, ∇)
    zₖₗ = directional_derivative(tri, k, ℓ, N₀, ∇)
    zₖₘ = directional_derivative(tri, k, m, N₀, ∇)
    zₗᵢ = directional_derivative(tri, ℓ, i, N₀, ∇)
    zₗⱼ = directional_derivative(tri, ℓ, j, N₀, ∇)
    zₗₖ = directional_derivative(tri, ℓ, k, N₀, ∇)
    zₗₘ = directional_derivative(tri, ℓ, m, N₀, ∇)
    zₘᵢ = directional_derivative(tri, m, i, N₀, ∇)
    zₘⱼ = directional_derivative(tri, m, j, N₀, ∇)
    zₘₖ = directional_derivative(tri, m, k, N₀, ∇)
    zₘₗ = directional_derivative(tri, m, ℓ, N₀, ∇)
    zᵢⱼₖ = hessian_form(tri, i, j, k, N₀, H)
    zᵢⱼₗ = hessian_form(tri, i, j, ℓ, N₀, H)
    zᵢⱼₘ = hessian_form(tri, i, j, m, N₀, H)
    zᵢₖₗ = hessian_form(tri, i, k, ℓ, N₀, H)
    zᵢₖₘ = hessian_form(tri, i, k, m, N₀, H)
    zᵢₗₘ = hessian_form(tri, i, ℓ, m, N₀, H)
    zⱼᵢₗ = hessian_form(tri, j, i, ℓ, N₀, H)
    zⱼᵢₖ = hessian_form(tri, j, i, k, N₀, H)
    zᵢᵢₘ = hessian_form(tri, i, i, m, N₀, H)
    zⱼₖₗ = hessian_form(tri, j, k, ℓ, N₀, H)
    zⱼₖₘ = hessian_form(tri, j, k, m, N₀, H)
    zⱼₗₘ = hessian_form(tri, j, ℓ, m, N₀, H)
    zₖᵢⱼ = hessian_form(tri, k, i, j, N₀, H)
    zₖᵢₗ = hessian_form(tri, k, i, ℓ, N₀, H)
    zₖᵢₘ = hessian_form(tri, k, i, m, N₀, H)
    zₖⱼₗ = hessian_form(tri, k, j, ℓ, N₀, H)
    zₖⱼₘ = hessian_form(tri, k, j, m, N₀, H)
    zₖₗₘ = hessian_form(tri, k, ℓ, m, N₀, H)
    zₗᵢⱼ = hessian_form(tri, ℓ, i, j, N₀, H)
    zₗᵢₖ = hessian_form(tri, ℓ, i, k, N₀, H)
    zₗᵢₘ = hessian_form(tri, ℓ, i, m, N₀, H)
    zₗⱼₖ = hessian_form(tri, ℓ, j, k, N₀, H)
    zₗⱼₘ = hessian_form(tri, ℓ, j, m, N₀, H)
    zₗₖₘ = hessian_form(tri, ℓ, k, m, N₀, H)
    zₘᵢⱼ = hessian_form(tri, m, i, j, N₀, H)
    zₘᵢₖ = hessian_form(tri, m, i, k, N₀, H)
    zₘᵢₗ = hessian_form(tri, m, i, ℓ, N₀, H)
    zₘⱼₖ = hessian_form(tri, m, j, k, N₀, H)
    zₘⱼₗ = hessian_form(tri, m, j, ℓ, N₀, H)
    zₘₖₗ = hessian_form(tri, m, k, ℓ, N₀, H)
    fᵢⱼₖₗₘ = @. (zᵢ + zⱼ + zₖ + zₗ + zₘ) / 5 +
                (zᵢⱼ + zᵢₖ + zᵢₗ + zᵢₘ + zⱼᵢ + zⱼₖ + zⱼₗ + zⱼₘ + zₖᵢ + zₖⱼ + zₖₗ +
                 zₖₘ + zₗᵢ + zₗⱼ + zₗₖ + zₗₘ + zₘᵢ + zₘⱼ + zₘₖ + zₘₗ) / 30 +
                (zᵢⱼₖ + zᵢⱼₗ + zᵢⱼₘ + zᵢₖₗ + zᵢₖₘ + zᵢₗₘ + zⱼᵢₗ + zⱼᵢₖ + zᵢᵢₘ +
                 zⱼₖₗ + zⱼₖₘ + zⱼₗₘ + zₖᵢⱼ + zₖᵢₗ + zₖᵢₘ + zₖⱼₗ + zₖⱼₘ + zₖₗₘ +
                 zₗᵢⱼ + zₗᵢₖ + zₗᵢₘ + zₗⱼₖ + zₗⱼₘ + zₗₖₘ + zₘᵢⱼ + zₘᵢₖ + zₘᵢₗ +
                 zₘⱼₖ + zₘⱼₗ + zₘₖₗ) / 180
    return fᵢⱼₖₗₘ
end

function get_contrib(tri, i, j, k, ℓ, m, N₀, ∇, H, z, case)
    if case == 1
        fᵢᵢᵢᵢᵢ = _hiyoshi_case_1(m, N₀, z)
        return fᵢᵢᵢᵢᵢ, 120
    elseif case == 2
        fᵢᵢᵢᵢⱼ = _hiyoshi_case_2(tri, ℓ, m, N₀, ∇, z)
        return fᵢᵢᵢᵢⱼ, 24
    elseif case == 3
        fᵢᵢⱼⱼⱼ = _hiyoshi_case_3(tri, k, ℓ, N₀, ∇, H, z)
        return fᵢᵢⱼⱼⱼ, 12
    elseif case == 4
        fᵢᵢᵢⱼₖ = _hiyoshi_case_4(tri, k, ℓ, m, N₀, ∇, H, z)
        return fᵢᵢᵢⱼₖ, 6
    elseif case == 5
        fᵢᵢⱼⱼₖ = _hiyoshi_case_5(tri, i, k, m, N₀, ∇, H, z)
        return fᵢᵢⱼⱼₖ, 4
    elseif case == 6
        fᵢᵢⱼₖₗ = _hiyoshi_case_6(tri, j, k, ℓ, m, N₀, ∇, H, z)
        return fᵢᵢⱼₖₗ, 2
    else
        fᵢⱼₖₗₘ = _hiyoshi_case_7(tri, i, j, k, ℓ, m, N₀, ∇, H, z)
        return fᵢⱼₖₗₘ, 1
    end
end

function compute_natural_coordinates(::Hiyoshi, tri, interpolation_point, cache=NaturalNeighboursCache(tri); kwargs...)
    return _compute_sibson_coordinates(tri, interpolation_point, cache; kwargs...)
end