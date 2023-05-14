struct NaturalCoordinates{F,I,T<:Triangulation}
    coordinates::Vector{F}
    indices::Vector{I}
    interpolation_point::NTuple{2,F}
    triangulation::T
end
function Base.show(io::IO, ::MIME"text/plain", nc::NaturalCoordinates{F,I}) where {F,I}
    coordinates = get_coordinates(nc)
    indices = get_indices(nc)
    interpolation_point = get_interpolation_point(nc)
    println(io, "NaturalCoordinates{", F, ",", I, "}")
    println(io, "    u: ", interpolation_point)
    println(io, "    λ: ", coordinates)
    print(io, "    k: ", indices)
end
get_coordinates(nc::NaturalCoordinates) = nc.coordinates
get_indices(nc::NaturalCoordinates) = nc.indices
get_interpolation_point(nc::NaturalCoordinates) = nc.interpolation_point
get_triangulation(nc::NaturalCoordinates) = nc.triangulation

