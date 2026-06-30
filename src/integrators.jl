
abstract type Integrator end
struct Euler <: Integrator end
struct Leapfrog <: Integrator end
struct RK4 <: Integrator end
struct VelocityVerlet <: Integrator end

# Euler integrator
function step!(integrator::Euler, system::System, dt::Float64)
  accels = get_accelerations(system)

  for (body, a) ∈ zip(system.bodies, accels)
    old_v = body.velocity
    body.velocity += a * dt
    body.position += old_v * dt
  end
end

# Leapfrog
function step!(integrator::Leapfrog, system::System, dt::Float64)
  dt2 = dt / 2
  accels = get_accelerations(system)
  for (body, a) ∈ zip(system.bodies, accels)
    body.velocity += a * dt2
    body.position += body.velocity * dt
  end

  accels = get_accelerations(system)
  for (body, a) ∈ zip(system.bodies, accels)
    body.velocity += a * dt2
  end
end

# RK4
function step!(integrator::RK4, system::System, dt::Float64)
  dt2 = dt / 2
  dt6 = dt / 6
  accels = get_accelerations(system)
  k1_v = SVector{3,Float64}[]
  k1_x = SVector{3,Float64}[]
  system2 = System(Body[])
  for (body, a) ∈ zip(system.bodies, accels)
    push!(k1_v, a)
    push!(k1_x, body.velocity)

    new_body = Body(body.position + body.velocity * dt2, body.velocity, body.mass)
    push!(system2.bodies, new_body)
  end

  accels = get_accelerations(system2)
  k2_v = SVector{3,Float64}[]
  k2_x = SVector{3,Float64}[]
  system3 = System(Body[])
  for (body, a, v) ∈ zip(system.bodies, accels, k1_v)
    new_vel = body.velocity + dt2 * v
    push!(k2_v, a)
    push!(k2_x, new_vel)

    new_body = Body(body.position + dt2 * new_vel, new_vel, body.mass)
    push!(system3.bodies, new_body)
  end

  accels = get_accelerations(system3)
  k3_v = SVector{3,Float64}[]
  k3_x = SVector{3,Float64}[]
  system4 = System(Body[])
  for (body, a, v) ∈ zip(system.bodies, accels, k2_v)
    new_vel = body.velocity + dt2 * v
    push!(k3_v, a)
    push!(k3_x, new_vel)

    new_body = Body(body.position + dt2 * new_vel, new_vel, body.mass)
    push!(system4.bodies, new_body)
  end

  accels = get_accelerations(system4)
  k4_v = SVector{3,Float64}[]
  k4_x = SVector{3,Float64}[]
  for (body, a, v) ∈ zip(system.bodies, accels, k3_v)
    push!(k4_v, a)
    push!(k4_x, body.velocity + dt * v)
  end

  i = 1
  for body ∈ system.bodies
    body.velocity += dt6 * (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i])
    body.position += dt6 * (k1_x[i] + 2*k2_x[i] + 2*k3_x[i] + k4_x[i])
    i += 1
  end
end

# Velocity Verlet
function step!(integrator::VelocityVerlet, system::System, dt::Float64)
  accels = get_accelerations(system)

  for (body, a) ∈ zip(system.bodies, accels)
    body.position += body.velocity * dt + 0.5 * a * (dt^2)
  end

  accels_next = get_accelerations(system)
  for (body, a, a_next) ∈ zip(system.bodies, accels, accels_next)
    body.velocity += 0.5 * (a + a_next) * dt
  end
end