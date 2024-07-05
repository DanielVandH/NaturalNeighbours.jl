using NaturalNeighbours
using Test
using SafeTestsets

include("jet_aqua.jl")

@testset "Interpolation" begin
    @testset "Coordinates" begin
        @info "Testing NaturalCoordinates"
        @safetestset "Natural Coordinates" begin
            include("interpolation/coordinates/natural_coordinates.jl")
        end
        @info "Testing Utils"
        @safetestset "Utils" begin
            include("interpolation/coordinates/utils.jl")
        end
    end
    @info "Testing Basic Tests"
    @safetestset "Basic Tests" begin
        include("interpolation/basic_tests.jl")
    end
    @info "Testing Precision"
    @safetestset "Precision" begin
        include("interpolation/precision.jl")
    end
    @info "Testing Extrapolation"
    @safetestset "Extrapolation" begin
        include("interpolation/extrapolation.jl")
    end
    @info "Testing Structs"
    @safetestset "Structs" begin
        include("interpolation/structs.jl")
    end
    @info "Testing Influence"
    @safetestset "Influence" begin
        include("interpolation/influence.jl")
    end 
    @info "Testing Constrained Triangulation"
    @safetestset "Constrained Triangulations" begin
        include("interpolation/constrained.jl")
    end
end

@testset "Differentiation" begin
    @info "Testing Basic Tests"
    @safetestset "Basic Tests" begin
        include("differentiation/basic_tests.jl")
    end
    @info "Testing Structs"
    @safetestset "Structs" begin
        include("differentiation/structs.jl")
    end
    @info "Testing Utils"
    @safetestset "Utils" begin
        include("differentiation/utils.jl")
    end
end

@testset "Documentation Examples" begin
    @info "Testing README Example"
    @safetestset "README Example" begin
        include("doc_examples/readme_example.jl")
    end
    @info "Testing Interpolation Example"
    @safetestset "Interpolation Example" begin
        include("doc_examples/interpolation.jl")
    end
    @info "Testing Differentiation Example"
    @safetestset "Differentiation Example" begin
        include("doc_examples/differentiation.jl")
    end
    @info "Testing Interpolation Math"
    @safetestset "Interpolation Math" begin
        include("doc_examples/interpolation_math.jl")
    end
    @info "Testing Switzerland"
    @safetestset "Switzerland" begin
        include("doc_examples/swiss.jl")
    end
    #if get(ENV, "CI", "false") == "false"
    #    @safetestset "Comparison" begin
    #        include("doc_examples/interpolant_comparisons.jl")
    #    end
    #end
end

using NaturalNeighbours
using PyPlot
using BenchmarkTools

function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]
    X, Y
end

f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
const x = vec([(i - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
const y = vec([(j - 1) / 9 for i in (1, 3, 4, 5, 8, 9, 10), j in (1, 2, 3, 5, 6, 7, 9, 10)])
const z = f.(x, y)

const xg = LinRange(0, 1, 100)
const yg = LinRange(0, 1, 100)

function nn_benchmark()
    itp = interpolate(x, y, z; derivatives=false)
    X, Y = meshgrid(xg, yg)
    _x, _y = vec(X), vec(Y)
    return itp(_x, _y; method=Triangle(; allow_cache = true), parallel=false)
end
function mpl_benchmark()
    triang = matplotlib.tri.Triangulation(x, y)
    interpolator = matplotlib.tri.LinearTriInterpolator(triang, z)
    X, Y = meshgrid(xg, yg)
    return interpolator(X, Y)
end

@benchmark nn_benchmark()
#=
julia> @benchmark nn_benchmark()
BenchmarkTools.Trial: 297 samples with 1 evaluation.
 Range (min … max):  15.082 ms … 95.031 ms  ┊ GC (min … max): 0.00% … 46.52%
 Time  (median):     15.459 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   16.848 ms ±  7.385 ms  ┊ GC (mean ± σ):  4.69% ±  7.69%

  █▄  
  ██▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▅▄▄▄▅ ▅
  15.1 ms      Histogram: log(frequency) by time      48.8 ms <

 Memory estimate: 1.85 MiB, allocs estimate: 48720.
=#
@benchmark mpl_benchmark()
#=
julia> @benchmark mpl_benchmark()
BenchmarkTools.Trial: 2832 samples with 1 evaluation.
 Range (min … max):  1.030 ms …  14.261 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.822 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.759 ms ± 417.625 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

       ▄                       ▄▅▆█▆▅▇▄▅▄▂▁▂▁▁▁
  ▂▄▅▆███▇▅▆▄▃▃▃▄▃▂▄▃▄▅▅▄▅▄▅▆▆█████████████████▆▆▅▅▄▄▃▃▂▂▂▂▂▂ ▄
  1.03 ms         Histogram: frequency by time        2.46 ms <

 Memory estimate: 237.89 KiB, allocs estimate: 74.
=#

# Your observed 16x difference is real. What if we use threads?

function nn_benchmark_threads() # I'm using 10 threads
    itp = interpolate(x, y, z; derivatives=false)
    X, Y = meshgrid(xg, yg)
    _x, _y = vec(X), vec(Y)
    return itp(_x, _y; method=Triangle(; allow_cache = true), parallel=true)
end

@benchmark nn_benchmark_threads()

#=
julia> @benchmark nn_benchmark_threads()
BenchmarkTools.Trial: 1322 samples with 1 evaluation.
 Range (min … max):  2.352 ms … 94.722 ms  ┊ GC (min … max):  0.00% … 80.89%
 Time  (median):     2.645 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   3.770 ms ±  7.382 ms  ┊ GC (mean ± σ):  19.06% ±  9.26%

  █▁ 
  ██▅▁▄▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▆▇ ▇
  2.35 ms      Histogram: log(frequency) by time     43.7 ms <

 Memory estimate: 1.84 MiB, allocs estimate: 47709.
=#






