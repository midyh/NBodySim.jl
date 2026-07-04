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

@testset "Diagnostic functions" begin
  test_system = make_test_system()
  momentum = total_momentum(test_system)
  angular = total_angular_momentum(test_system)
  energy = total_energy(test_system)

  momentum_hand = SA[0.0, 2*(1/sqrt(3))+1*(-2/sqrt(3)), 0.0]
  @test momentum == momentum_hand

  angular_hand = #=2 * SA[0.0, 0.0, 0.0] + =#1 * SA[0.0, 0.0, -2/sqrt(3)]
  @test angular == angular_hand

  energy_hand = 0.5 * 2 * (1/sqrt(3))^2 + 0.5 * 1 * (2/sqrt(3))^2 - (1 * 2 * 1)/1
  @test energy ≈ energy_hand
end

@testset "Euler Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  initial_energy = total_energy(test_system)
  println(initial_angular)
  momentums, angulars, energies = simulate_with_diagnostics(test_system, Euler(), 0.001, 100)

  for m ∈ momentums
    @test initial_momentum ≈ m
  end

  # Euler has a monotonically increasing drift, growing by roughly 3.46 * 10^(-4)
  # Proportional to dt. Here we verify that the drift is consistent
  @test norm(initial_angular - angulars[end]) > 1e-6

  @test abs(energies[end] - initial_energy) > 1e-4
end

@testset "Leapfrog Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  initial_energy = total_energy(test_system)
  momentums, angulars, energies = simulate_with_diagnostics(test_system, Leapfrog(), 0.01, 100)

  for (m, a) ∈ zip(momentums, angulars)
    @test initial_momentum ≈ m
    @test initial_angular ≈ a
  end

  @test all(abs.(E .- initial_energy) .< 1e-4 for E in energies)
end

@testset "Velocity Verlet Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  initial_energy = total_energy(test_system)
  momentums, angulars, energies = simulate_with_diagnostics(test_system, VelocityVerlet(), 0.01, 100)

  for (m, a) ∈ zip(momentums, angulars)
    @test initial_momentum ≈ m
    @test initial_angular ≈ a
  end

  @test all(abs.(E .- initial_energy) .< 1e-4 for E in energies)
end

@testset "RK4 Integrator" begin
  test_system = make_test_system()

  initial_momentum = total_momentum(test_system)
  initial_angular = total_angular_momentum(test_system)
  initial_energy = total_energy(test_system)
  momentums, angulars, energies = simulate_with_diagnostics(test_system, RK4(), 0.01, 100)

  for m ∈ momentums
    @test initial_momentum ≈ m
  end

  # Similar to Euler, however rate of drift is 12 times slower than Euler
  @test norm(initial_angular - angulars[end]) > 1e-3

  @test abs(energies[end] - initial_energy) > 1e-3
end
