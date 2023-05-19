struct NaturalNeighboursDifferentiator{I,O}
    interpolant::I 
    function NaturalNeighboursDifferentiator(itp::I, order) where {I}
        @assert order ∈ (1, 2) "order must be 1 or 2."
        return new{I, order}(itp)
    end
end
get_interpolant(d::NaturalNeighboursDifferentiator) = d.interpolant
function Base.show(io::IO, ::MIME"text/plain", d::NaturalNeighboursDifferentiator{I,O}) where {I,O}
    z = get_z(get_interpolant(d))
    ∇ = get_gradient(get_interpolant(d))
    ℋ = get_hessian(get_interpolant(d))
    println(io, "Natural Neighbour Differentiator")
    println(io, "    Order: ", O)
    println(io, "    z: ", z)
    println(io, "    ∇: ", ∇)
    print(io, "    H: ", ℋ)
end