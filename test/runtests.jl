using Test
using LinearAlgebra
using NBodySim

include("../src/includes.jl")

make_body1() = Body(
  SA[0.0, 0.0, 0.0],
  SA[0.0, 1/sqrt(3), 0.0],
  2.0,
)

make_body2() = Body(
  SA[1.0, 0.0, 0.0],
  SA[0.0, -2/sqrt(3), 0.0],
  1.0,
)

make_test_system() = System([
  make_body1(),
  make_body2(),
])

function simulate_with_diagnostics(system::System, integrator::Integrator, dt::Float64, n_steps::Int64)
  momentums = SVector{3,Float64}[]
  angulars = SVector{3,Float64}[]

  for _ ∈ 1:n_steps
    step!(integrator, system, dt)
    push!(momentums, total_momentum(system))
    push!(angulars, total_angular_momentum(system))
  end

  return momentums, angulars
end

@testset "Euler Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  println(initial_angular)
  momentums, angulars = simulate_with_diagnostics(test_system, Euler(), 0.001, 100)

  for m ∈ momentums
    @test initial_momentum ≈ m
  end

  # Euler has a monotonically increasing drift, growing by roughly 3.46 * 10^(-4)
  # Proportional to dt. Here we verify that the drift is consistent
  @test norm(initial_angular - angulars[end]) > 1e-6
end

@testset "Leapfrog Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  momentums, angulars = simulate_with_diagnostics(test_system, Leapfrog(), 0.01, 100)

  for (m, a) ∈ zip(momentums, angulars)
    @test initial_momentum ≈ m
    @test initial_angular ≈ a
  end
end

@testset "Velocity Verlet Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  momentums, angulars = simulate_with_diagnostics(test_system, VelocityVerlet(), 0.01, 100)

  for (m, a) ∈ zip(momentums, angulars)
    @test initial_momentum ≈ m
    @test initial_angular ≈ a
  end
end

@testset "RK4 Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  momentums, angulars = simulate_with_diagnostics(test_system, RK4(), 0.01, 100)

  for m ∈ momentums
    @test initial_momentum ≈ m
  end

  # Similar to Euler, however rate of drift is 12 times slower than Euler
  @test norm(initial_angular - angulars[end]) > 1e-3
end
