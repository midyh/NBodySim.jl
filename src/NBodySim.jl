module NBodySim

using StaticArrays

export Body, System, SA
export simulate, animate
export Euler, Leapfrog, VelocityVerlet, RK4Z

include("includes.jl")

end