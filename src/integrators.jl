
abstract type Integrator end
struct Euler <: Integrator end

function step!(integrator::Euler, system::System, dt::Float64)
  accels = get_accelerations(system)

  for (body::Body, a::SVector) ∈ zip(system.bodies, accels)
    old_v = body.velocity
    body.velocity += a * dt
    body.position += old_v * dt
  end
end