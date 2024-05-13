function get_taylor_neighbourhood!(S, S′, tri, i::Integer, d, ::Any; rng=nothing)
    iterated_neighbourhood!(S, tri, i, d)
    return one(i), S
end

function get_taylor_neighbourhood!(S, S′, tri, p, d, n_cache=NaturalNeighboursCache(tri); rng=Random.default_rng())
    @assert d ∈ (1, 2, 3) "d must be 1, 2 or 3."
    nc = compute_natural_coordinates(Sibson(), tri, p, n_cache; rng)
    return get_taylor_neighbourhood!(S, S′, tri, d, nc)
end

function get_taylor_neighbourhood!(S, S′, tri, d, nc::NaturalCoordinates)
    coordinates = get_coordinates(nc)
    envelope = get_indices(nc)
    idx = findfirst(isone, coordinates)
    if !isnothing(idx)
        return get_taylor_neighbourhood!(S, S′, tri, envelope[idx], d, nc)
    elseif d == 2
        # Bit hard to use iterated_neighbourhood here, since we want to make sure we don't overwrite any 
        # coordinates from the initial set of natural neighbours. So we just do it manually. Will
        # need to generalise if we ever want higher order derivatives.
        empty!(S′)
        for i in envelope
            for j in get_neighbours(tri, i)
                if !is_ghost_vertex(j) && j ∉ envelope
                    push!(S′, j)
                end
            end
        end
        for s in S′
            push!(envelope, s)
        end
    elseif d == 3
        empty!(S′)
        for i in envelope
            for j in get_neighbours(tri, i)
                if !is_ghost_vertex(j) && j ∉ envelope
                    push!(S′, j)
                end
            end
        end
        empty!(S)
        for s in S′
            for j in get_neighbours(tri, s)
                if !is_ghost_vertex(j) && j ∉ envelope && j ∉ S′
                    push!(S, j)
                end
            end
        end
        for s in S
            push!(envelope, s)
        end
        for s in S′
            push!(envelope, s)
        end
    end
    return coordinates, envelope
end

function get_λ(coordinates::Vector{F}, j, use_sibson_weight) where {F}
    if !use_sibson_weight || j > lastindex(coordinates)
        return one(F)
    else
        return coordinates[j]
    end
end
get_λ(coordinates::F, j, use_sibson_weight) where {F<:Number} = one(F)

setval!(A::Vector, i, val) = A[i] = val
setval!(A::Matrix, i, val) = A[:, i] .= val