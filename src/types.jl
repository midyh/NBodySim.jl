import StaticArrays

mutable struct Body
  position::StaticArrays.SVector{3,Float64}
  velocity::StaticArrays.SVector{3,Float64}
  mass::Float64
end

struct System
  bodies::Vector{Body}
end
