import DelaunayTriangulation: DelaunayTriangulation

import ChunkSplitters: chunks 
import ElasticArrays 
import LinearAlgebra
import Random 

#### Interpolators 
"""
    AbstractInterpolator

Abstract type for defining interpolants. Within this package we define 
the following interpolants:

- [`Sibson0`](@ref)
- [`Sibson1`](@ref)
- [`Triangle`](@ref)
- [`Nearest`](@ref)
- [`Laplace`](@ref)
- [`Farin`](@ref)
- [`Hiyoshi2`](@ref)
"""
abstract type AbstractInterpolator end 

"""
    Sibson0() <: AbstractInterpolator 

Defines the Sibson interpolant with `C(0)` continuity at the data sites. 
"""
struct Sibson0 <: AbstractInterpolator end 

"""
    Sibson1() <: AbstractInterpolator

Defines the Sibson interpolant with `C(1)` continuity at the data sites.
"""
struct Sibson1 <: AbstractInterpolator end 

"""
    Triangle() <: AbstractInterpolator

Defines the triangle interpolant, an interpolant with `C(0)` continity 
at the data sites that is piecewise linear within each triangle of 
the associated triangulation.
"""
struct Triangle <: AbstractInterpolator end 

"""
    Nearest() <: AbstractInterpolator

Defines the nearest neighbour interpolant, an interpolant which maps each point 
to the value of the function at the point's nearest data site.
"""
struct Nearest <: AbstractInterpolator end 

"""
    Laplace() <: AbstractInterpolator

Deines Laplace's interpolant with `C(0)` continity at the data sites.
"""
struct Laplace <: AbstractInterpolator end 

"""
    Farin() <: AbstractInterpolator

Defines Farin's interpolant with `C(1)` continity at the data sites.
"""
struct Farin <: AbstractInterpolator end 

"""
    Hiyoshi2() <: AbstractInterpolator

Defines Hiyoshi's interpolant with `C(2)` continuity at the data sites.
"""
struct Hiyoshi2 <: AbstractInterpolator end 

#### Utilities 
@inline _zero(::Type{F}) where {F} = zero(F)
@inline _zero(::Type{NTuple{N,F}}) where {N,F} = ntuple(i -> _zero(F), Val(N))
@inline zero!(A, i) = setval!(A, i, _zero(eltype(A)))

@inline _copy(z) = copy(z)
@inline _copy(z::Tuple) = z

#### Data interface 

@doc """
    ndata(z) -> Int 

Returns the number of data points associated with the data array `z`. For vectors,
returns `length(z)`, and for matrices returns `size(z, 2)`.
"""
ndata
@inline ndata(z::AbstractVector) = length(z)
@inline ndata(z::AbstractMatrix) = size(z, 2)

@doc """
    fdim(z) -> Int 

Returns the dimension of the function that generated the data `z`. For 
vectors of numbers, returns `1`, and for matrices of numbers returns `size(z, 1)`.
"""
fdim
@inline fdim(z::AbstractVector{<:Number}) = 1
@inline fdim(z::AbstractVector{<:AbstractVector}) = length(get_data(z, 1))
@inline fdim(z::AbstractMatrix{<:Number}) = size(z, 1)

@doc """
    zrange(z)

Returns an iterator over the indices of the data points `z`. For vectors of numbers, 
returns `eachindex(z)`, and for matrices of numbers returns `axes(z, 2)`.
"""
zrange
@inline zrange(z::AbstractVector{<:Number}) = eachindex(z)
@inline zrange(z::AbstractVector{<:AbstractVector}) = eachindex(get_data(z, 1))
@inline zrange(z::AbstractMatrix{<:Number}) = axes(z, 2)

@doc """
    setval!(A, i, val)

Sets the `i`th data point in `A` to be `val`. This is essentially like `setindex!`, except 
we are thinking in terms of data points rather than indices. For instance, for vectors of numbers 
this definition is `A[i] = val`, while for matrices of numbers we use `A[:, i] .= val`.
"""
setval!
@inline setval!(A::AbstractVector{<:Number}, i, val) = A[i] = val 
@inline setval!(A::AbstractVector{<:AbstractVector}, i, val) = A[i] .= val
@inline setval!(A::AbstractMatrix{<:Number}, i, val) = A[:, i] .= val 

@doc """
    get_data(z, i)

Returns the `i`th data point from `z`. For vectors, returns `z[i]`,
and for matrices it returns `view(z, :, i)`.
"""
get_data
@inline get_data(z::AbstractVector, i) = z[i]
@inline get_data(z::AbstractMatrix{<:Number}, i) = view(z, :, i)

"""
    is_scalar(z) -> Bool 

Tests if the data points associated with `z` are scalar. This is tested 
using [`fdim`], i.e. `is_scalar(z)` is the same as `fdim(z) == 1`.
"""
is_scalar
@inline is_scalar(z::AbstractVector{<:Number}) = true 
@inline is_scalar(z::Union{AbstractVector{<:AbstractVector}, AbstractMatrix{<:Number}}) = false 
@inline is_scalar(z) = fdim(z) == 1 

#### Caches 
struct NaturalNeighboursCache{F,I,H,E,R}
    coordinates::Vector{F}
    envelope::Vector{I}
    insertion_event_history::H
    poly_points::Vector{NTuple{2,F}}
    temp_adjacent::Adjacent{I,E}
    last_triangle::R
end
function NaturalNeighboursCache(tri::Triangulation)
    F, I, E = number_type(tri), integer_type(tri), edge_type(tri)
    coordinates = F[]
    envelope = I[]
    insertion_event_history = InsertionEventHistory(tri)
    poly_points = NTuple{2,number_type(tri)}[]
    temp_adjacent = Adjacent{I,E}()
    last_triangle = Ref(triangle_vertices(first(each_solid_triangle(tri))))
    return NaturalNeighboursCache(coordinates, envelope, insertion_event_history, poly_points, temp_adjacent, last_triangle)
end

struct DerivativeCache{I, F, T}
    iterated_neighbourhood::Set{I}
    second_iterated_neighbourhood::Set{I}
    linear_matrix::ElasticMatrix{F, Vector{F}}
    quadratic_matrix::ElasticMatrix{F, Vector{F}}
    quadratic_matrix_no_cubic::ElasticMatrix{F, Vector{F}}
    rhs_vector::T 
    gradient::Matrix{F} 
    hessian::Matrix{F}
end
function DerivativeCache(tri::Triangulation, z)
    isscalar = is_scalar(z)
    I, F = integer_type(tri), number_type(tri) 
    iterated_neighbourhood = Set{I}()
    second_iterated_neighbourhood = Set{I}()
    linear_matrix = ElasticMatrix{F, Vector{F}}(undef, 2, 0)
    quadratic_matrix = ElasticMatrix{F, Vector{F}}(undef, 9, 0)
    quadratic_matrix_no_cubic = ElasticMatrix{F, Vector{F}}(undef, 5, 0)
    rhs_vector = isscalar ? zeros(F, 0) : ElasticMatrix{F, Vector{F}}(undef, fdim(z), 0)
    gradient = Matrix{F}(undef, 2, fdim(z))
    hessian = Matrix{F}(undef, 3, fdim(z))
    sizehint_m = 2^5
    sizehint!(iterated_neighbourhood, sizehint_m)
    sizehint!(second_iterated_neighbourhood, sizehint_m)
    sizehint!(linear_matrix, (2, sizehint_m))
    sizehint!(quadratic_matrix, (9, sizehint_m))
    sizehint!(quadratic_matrix_no_cubic, (5, sizehint_m))
    isscalar ? sizehint!(rhs_vector, sizehint_m) : sizehint!(rhs_vector, (fdim(z), sizehint_m))
    return DerivativeCache(iterated_neighbourhood, second_iterated_neighbourhood, linear_matrix, quadratic_matrix, quadratic_matrix_no_cubic, rhs_vector, gradient, hessian)
end

