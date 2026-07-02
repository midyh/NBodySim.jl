using LinearAlgebra

function total_momentum(system::System)::SVector{3,Float64}
  p_total = zeros(SVector{3,Float64})
  for body ∈ system.bodies
    p_total += body.mass * body.velocity
  end

  return p_total
end

function total_angular_momentum(system::System)::SVector{3,Float64}
  L_total = zeros(SVector{3,Float64})

  for body ∈ system.bodies
    L_total += body.mass * cross(body.position, body.velocity)
  end

  return L_total
end

function total_energy(system::System, G::Float64=1.0)::Float64
  E_total = 0
  for body ∈ system.bodies
    E_total += 0.5 * body.mass * norm(body.velocity)^2
    for body2 ∈ system.bodies
      if body === body2
        continue
      end

      E_total += (-G * body.mass * body2.mass) / (norm(body.position - body2.position))
    end
  end

  return E_total / 2
end