using CairoMakie
using BenchmarkTools

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

function plot_energies(n_steps)
  _, _, euler = simulate_with_diagnostics(make_test_system(), Euler(), 0.01, n_steps)
  _, _, leapfrog = simulate_with_diagnostics(make_test_system(), Leapfrog(), 0.01, n_steps)
  _, _, velover = simulate_with_diagnostics(make_test_system(), VelocityVerlet(), 0.01, n_steps)
  _, _, rk4 = simulate_with_diagnostics(make_test_system(), RK4(), 0.01, n_steps)

  steps = 1:n_steps

  f = Figure()
  ax = Axis(f[1, 1],
    title="Energy drift over iterations (symplectic only)",
    xlabel="Iterations",
    ylabel="Energy",
  )
  lines!(ax, steps, leapfrog .- leapfrog[1], label="Leapfrog", color=:red)
  lines!(ax, steps, velover .- velover[1], label="Velocity Verlet", color=:green)
  axislegend(ax, position=:lt)

  f1 = Figure()
  ax1 = Axis(f1[1, 1],
    title="Energy drift over iterations",
    xlabel="Iterations",
    ylabel="Energy",
  )
  lines!(ax1, steps, euler .- euler[1], label="Euler")
  lines!(ax1, steps, leapfrog .- leapfrog[1], label="Leapfrog", color=:red)
  lines!(ax1, steps, rk4 .- rk4[1], label="RK4")
  axislegend(ax1, position=:lt)

  display(f)
  display(f1)
end

plot_energies(1000)

function bench_performance(n_steps)
  euler_system = make_test_system()
  leapfrog_system = make_test_system()
  velver_system = make_test_system()
  rk4_system = make_test_system()

  println("Euler")
  display(@benchmark simulate_with_diagnostics($euler_system, Euler(), 0.01, $n_steps))

  println("\nLeapfrog")
  display(@benchmark simulate_with_diagnostics($leapfrog_system, Leapfrog(), 0.01, $n_steps))

  println("\nVelocity Verlet")
  display(@benchmark simulate_with_diagnostics($velver_system, VelocityVerlet(), 0.01, $n_steps))

  println("\nRK4")
  display(@benchmark simulate_with_diagnostics($rk4_system, RK4(), 0.01, $n_steps))
end

results = bench_performance(100)