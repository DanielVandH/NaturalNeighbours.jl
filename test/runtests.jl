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