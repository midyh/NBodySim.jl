module NBodySim

using StaticArrays

include("includes.jl")

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

end