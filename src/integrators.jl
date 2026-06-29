
abstract type Integrator end
struct Euler <: Integrator end
struct Leapfrog <: Integrator end

function step!(integrator::Euler, system::System, dt::Float64)
  accels = get_accelerations(system)

  for (body, a) ∈ zip(system.bodies, accels)
    old_v = body.velocity
    body.velocity += a * dt
    body.position += old_v * dt
  end
end

function step!(integrator::Leapfrog, system::System, dt::Float64)
  accels = get_accelerations(system)
  for (body, a) ∈ zip(system.bodies, accels)
    body.velocity += a * (dt / 2)
    body.position += body.velocity * dt
  end

  accels = get_accelerations(system)
  for (body, a) ∈ zip(system.bodies, accels)
    body.velocity += a * (dt / 2)
  end
end