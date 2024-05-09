struct NaturalNeighboursCache{F,I,H,E,R}
    coordinates::Vector{F}
    envelope::Vector{I}
    insertion_event_history::H
    poly_points::Vector{NTuple{2,F}}
    temp_adjacent::Adjacent{I,E}
    last_triangle::R
end
get_coordinates(cache::NaturalNeighboursCache) = cache.coordinates
get_envelope(cache::NaturalNeighboursCache) = cache.envelope
get_insertion_event_history(cache::NaturalNeighboursCache) = cache.insertion_event_history
get_poly_points(cache::NaturalNeighboursCache) = cache.poly_points
get_temp_adjacent(cache::NaturalNeighboursCache) = cache.temp_adjacent
get_last_triangle(cache::NaturalNeighboursCache) = cache.last_triangle
function NaturalNeighboursCache(tri::Triangulation{P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE}) where {P,T,BN,W,I,E,Es,BC,BEM,GVM,GVR,BPL,C,BE}
    coordinates = number_type(tri)[]
    envelope = I[]
    insertion_event_history = InsertionEventHistory(tri)
    poly_points = NTuple{2,number_type(tri)}[]
    temp_adjacent = Adjacent{I,E}()
    last_triangle = (Ref ∘ triangle_vertices ∘ first ∘ each_solid_triangle)(tri)
    return NaturalNeighboursCache(coordinates, envelope, insertion_event_history, poly_points, temp_adjacent, last_triangle)
end