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

function simulate_with_diagnostics(system::System, integrator::Integrator, dt::Float64, n_steps::Int64)
  momentums = SVector{3,Float64}[]
  angulars = SVector{3,Float64}[]
  energies = Float64[]

  for _ ∈ 1:n_steps
    step!(integrator, system, dt)
    push!(momentums, total_momentum(system))
    push!(angulars, total_angular_momentum(system))
    push!(energies, total_energy(system))
  end

  return momentums, angulars, energies
end

function simulate_for_visualization(
  system::System,
  integrator::Integrator,
  dt::Float64,
  n_steps::Int64,
  G::Float64=1.0
)::Vector{Vector{SVector{3,Float64}}}
  all_positions = [SVector{3,Float64}[] for _ in system.bodies]

  for _ ∈ 1:n_steps
    step!(integrator, system, dt, G)

    for (i, body) in enumerate(system.bodies)
      push!(all_positions[i], body.position)
    end
  end

  return all_positions
end