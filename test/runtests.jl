using NaturalNeighbours
using Test
using SafeTestsets

@testset "Interpolation" begin
    @testset "Coordinates" begin
        @safetestset "Natural Coordinates" begin
            include("interpolation/coordinates/natural_coordinates.jl")
        end
        @safetestset "Utils" begin
            include("interpolation/coordinates/utils.jl")
        end
    end
    @safetestset "Basic Tests" begin
        include("interpolation/basic_tests.jl")
    end
    @safetestset "Precision" begin
        include("interpolation/precision.jl")
    end
    @safetestset "Extrapolation" begin
        include("interpolation/extrapolation.jl")
    end
    @safetestset "Structs" begin
        include("interpolation/structs.jl")
    end
end

@testset "Differentiation" begin
    @safetestset "Basic Tests" begin
        include("differentiation/basic_tests.jl")
    end
    @safetestset "Structs" begin
        include("differentiation/structs.jl")
    end
    @safetestset "Utils" begin
        include("differentiation/utils.jl")
    end
end

@testset "Documentation Examples" begin
    @safetestset "README Example" begin
        include("doc_examples/readme_example.jl")
    end
    @safetestset "Interpolation Example" begin
        include("doc_examples/interpolation.jl")
    end
    @safetestset "Differentiation Example" begin
        include("doc_examples/differentiation.jl")
    end
    @safetestset "Interpolation Math" begin
        include("doc_examples/interpolation_math.jl")
    end
    @safetestset "Switzerland" begin
        include("doc_examples/swiss.jl")
    end
end

@testset "Test Functions" begin
    @safetestset "Complete Analysis" begin
        include("test_functions/complete_analysis.jl")
    end
    @safetestset "Influence Regions" begin
        include("test_functions/influence.jl")
    end
end












