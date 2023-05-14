struct InterpolantCache{F,I,H,E,R}
    coordinates::Vector{F}
    envelope::Vector{I}
    insertion_event_history::H
    poly_points::Vector{NTuple{2,F}}
    temp_adjacent::Adjacent{I,E}
    last_triangle::R
end
get_coordinates(cache::InterpolantCache) = cache.coordinates
get_envelope(cache::InterpolantCache) = cache.envelope
get_insertion_event_history(cache::InterpolantCache) = cache.insertion_event_history
get_poly_points(cache::InterpolantCache) = cache.poly_points
get_temp_adjacent(cache::InterpolantCache) = cache.temp_adjacent
get_last_triangle(cache::InterpolantCache) = cache.last_triangle
function InterpolantCache(tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}
    coordinates = number_type(tri)[]
    envelope = I[]
    insertion_event_history = initialise_event_history(tri)
    poly_points = NTuple{2,number_type(tri)}[]
    temp_adjacent = Adjacent{I,E}()
    last_triangle = (Ref ∘ indices ∘ first ∘ each_solid_triangle)(tri)
    return InterpolantCache(coordinates, envelope, insertion_event_history, poly_points, temp_adjacent, last_triangle)
end