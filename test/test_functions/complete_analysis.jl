using ..NaturalNeighbours
const NNI = NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using LinearAlgebra
using ReferenceTests
using StatsBase
using DataStructures

to_mat(H::NTuple{3,Float64}) = [H[1] H[3]; H[3] H[2]]
to_mat(H) = H

include(normpath(@__DIR__, "../.", "helper_functions", "test_functions.jl"))

map(complete_test_function_analysis, [1,2,3,5,6])