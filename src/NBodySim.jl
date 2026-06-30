module NBodySim

using StaticArrays

include("types.jl")
include("dynamics.jl")
include("integrators.jl")

make_sun() = Body(
  SA[0.0, 0.0, 0.0],
  SA[0.0, 0.0, 0.0],
  10.0,
)

make_earth() = Body(
  SA[2.0, 2.0, 2.0],
  SA[0.0, 0.0, 0.0],
  2.0,
)

make_system() = System([
  make_sun(),
  make_earth(),
])

make_n_steps() = 5

make_dt() = 0.01

make_integrator() = Leapfrog()

main() = simulate(make_system(), make_integrator(), make_dt(), make_n_steps())

function simulate(system::System, integrator::Integrator, dt::Float64, n_steps::Int64)
  println("Using ", integrator)
  for i ∈ 1:n_steps
    step!(integrator, system, dt)
    println("Step ", i, ": ")
    for body ∈ system.bodies
      println(body.position)
    end
    println()
  end
end

end