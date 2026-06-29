using StaticArrays
using LinearAlgebra

function get_accelerations(system::System, G::Float64=1.0)
  accels = SVector{3,Float64}[]
  for b1::Body ∈ system.bodies
    acc = zeros(SVector{3,Float64})
    for b2::Body ∈ system.bodies
      if b1 === b2
        continue
      end

      dist_vec = b2.position - b1.position
      acc += (G * b2.mass * dist_vec) / norm(dist_vec)^3
    end
    push!(accels, acc)
  end
  return accels
end