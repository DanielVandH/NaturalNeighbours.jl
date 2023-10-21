using ..NaturalNeighbours
using DelaunayTriangulation
using CairoMakie
using ReferenceTests
using StableRNGs
using JET
using Aqua

Aqua.test_all(NaturalNeighbours; ambiguities=false, project_extras=false) # don't care about julia < 1.2
Aqua.test_ambiguities(NaturalNeighbours) # don't pick up Base and Core...

rng = StableRNG(123)
pts = [(rand(rng), rand(rng)) for _ in 1:50]
tri = triangulate(pts, rng=rng, delete_ghosts=false)
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
z = [f(x, y) for (x, y) in pts]
interpolate(tri, z; derivatives=false, parallel=false)

res = report_package(NaturalNeighbours; target_modules=(@__MODULE__,))
@test_opt target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=false, parallel=false)
@test_opt target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=false)
@test_opt target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=true)
@test_opt target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=true)
@test_opt target_modules = (@__MODULE__,) interpolate(tri, z)
itp = interpolate(tri, z)
@test_opt target_modules = (@__MODULE__,) differentiate(itp, 1)
@test_opt target_modules = (@__MODULE__,) differentiate(itp, 2)
∂ = differentiate(itp, 1)
@test_opt target_modules = (@__MODULE__,) ∂(0.3, 0.5; method = Iterative())
@test_opt target_modules = (@__MODULE__,) ∂(0.3, 0.5; method = Direct())
@test_opt target_modules = (@__MODULE__,) ∂(0.3, 0.5; method = Iterative(), project=false)
@test_opt target_modules = (@__MODULE__,) ∂(0.3, 0.5; method = Direct(),project=false)
@test_opt target_modules = (@__MODULE__,) ∂([0.3], [0.5], method=Iterative())
@test_opt target_modules = (@__MODULE__,) ∂([0.3], [0.5], method=Direct())
@test_opt target_modules = (@__MODULE__,) ∂([0.3], [0.5], method=Direct(), project=false)
@test_opt target_modules = (@__MODULE__,) ∂([0.3], [0.5], method=Iterative(),project=false)
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Sibson())
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Sibson(1))
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Farin())
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Hiyoshi())
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Laplace())
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Triangle())
@test_opt target_modules = (@__MODULE__,) itp(0.3, 0.5, method=Nearest())
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Sibson(), parallel=false)
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Sibson(1), parallel=true)
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Farin())
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Hiyoshi())
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Laplace())
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Triangle())
@test_opt target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Nearest())

@test_call target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=false, parallel=false)
@test_call target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=false)
@test_call target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=true)
@test_call target_modules = (@__MODULE__,) interpolate(tri, z; derivatives=true, parallel=true)
@test_call target_modules = (@__MODULE__,) interpolate(tri, z)
itp = interpolate(tri, z)
@test_call target_modules = (@__MODULE__,) differentiate(itp, 1)
@test_call target_modules = (@__MODULE__,) differentiate(itp, 2)
∂ = differentiate(itp, 1)
@test_call target_modules = (@__MODULE__,) ∂(0.3, 0.5)
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Sibson())
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Sibson(1))
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Farin())
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Hiyoshi())
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Laplace())
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Triangle())
@test_call target_modules = (@__MODULE__,) itp([0.3], [0.5], method=Nearest())