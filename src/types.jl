using StaticArrays

mutable struct Body
  position::SVector{3,Float64}
  velocity::SVector{3,Float64}
  mass::Float64
  Body(position, velocity, mass) = mass ≤ 0 ? error("mass cannot be 0 or less") : new(position, velocity, mass)
end

mutable struct System
  bodies::Vector{Body}
end
