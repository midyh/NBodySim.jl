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

  for i ∈ 1:length(system.bodies)
    E_total += 0.5 * system.bodies[i].mass * norm(system.bodies[i].velocity)^2
    for j ∈ (i+1):length(system.bodies)
      E_total += (-G * system.bodies[i].mass * system.bodies[j].mass) / (norm(system.bodies[i].position - system.bodies[j].position))
    end
  end

  return E_total
end