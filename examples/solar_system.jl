using CairoMakie

include("../src/includes.jl")

make_sun() = Body(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 0.0], 1.0)

make_mercury() = Body(
    SA[0.38709927, 0.0, 0.0], SA[0.0, 0.027650, 0.0], 1.651*1e-7
)

make_venus() = Body(
    SA[0.72333566, 0.0, 0.0], SA[0.0, 0.0202277, 0.0], 2.447*1e-6
)

make_earth() = Body(
    SA[1.00000261, 0.0, 0.0], SA[0.0, 0.017203, 0.0], 3.003*1e-6
)

make_mars() = Body(
    SA[1.52371034, 0.0, 0.0], SA[0.0, 0.013936, 0.0], 3.213 * 1e-7
)

make_solar_system() = System([make_sun(), make_mercury(), make_venus(), make_earth(), make_mars()])

function animate_solar_system(days::Int64=365)
    system = make_solar_system()
    pos = simulate(system, Leapfrog(), 1.0, days, (4*π^2)/(365.25)^2)

    fig = Figure()
    ax = Axis(fig[1, 1], title="Inner Solar System Orbit")
    hidespines!(ax, :t, :r)

    ax.xgridvisible = false
    ax.ygridvisible = false

    xs = Float64[]
    ys = Float64[]
    for traj ∈ pos, p ∈ traj
        push!(xs, p[1])
        push!(ys, p[2])
    end
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    limits!(ax, xmin - 0.05, xmax + 0.05, ymin - 0.05, ymax + 0.05)

    sun = Observable(Point2f0(0.0, 0.0))
    mercury = Observable(Point2f0(0.38709927, 0.0))
    mercury_orb = Observable(Point2f0[(0.38709927, 0.0)])
    venus = Observable(Point2f0(0.72333566, 0.0))
    venus_orb = Observable(Point2f0[(0.72333566, 0.0)])
    earth = Observable(Point2f0(1.00000261, 0.0))
    earth_orb = Observable(Point2f0[(1.00000261, 0.0)])
    mars = Observable(Point2f0(1.52371034, 0.0))
    mars_orb = Observable(Point2f0[(1.52371034, 0.0)])

    scatter!(ax, sun; color=:yellow, label="Sun", markersize=30, marker=:star5)
    scatter!(ax, mercury; color=:grey, label="Mercury")
    lines!(ax, mercury_orb; color=:grey)
    scatter!(ax, venus; color=:pink, label="Venus")
    lines!(ax, venus_orb; color=:pink)
    scatter!(ax, earth; color=:green, label="Earth")
    lines!(ax, earth_orb; color=:green)
    scatter!(ax, mars; color=:red, label="Mars")
    lines!(ax, mars_orb; color=:red)

    axislegend(ax, position=:lt)

    record(fig, "./docs/images/solar_system.mp4", 1:2:days; framerate=30) do step
        mercury_pos = Point2f0(pos[2][step][1], pos[2][step][2])
        venus_pos = Point2f0(pos[3][step][1], pos[3][step][2])
        earth_pos = Point2f0(pos[4][step][1], pos[4][step][2])
        mars_pos = Point2f0(pos[5][step][1], pos[5][step][2])

        mercury[] = mercury_pos
        venus[] = venus_pos
        earth[] = earth_pos
        mars[] = mars_pos

        mercury_orb[] = push!(mercury_orb[], mercury_pos)
        venus_orb[] = push!(venus_orb[], venus_pos)
        earth_orb[] = push!(earth_orb[], earth_pos)
        mars_orb[] = push!(mars_orb[], mars_pos)
    end
end

animate_solar_system()