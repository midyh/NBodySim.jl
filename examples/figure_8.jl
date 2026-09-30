using CairoMakie

include("../src/includes.jl")

make_body1() = Body(SA[-0.97000436, 0.24308753, 0.0], SA[0.93240737/2, 0.86473146/2, 0.0], 1.0)

make_body2() = Body(SA[0.97000436, -0.24308753, 0.0], SA[0.93240737/2, 0.86473146/2, 0.0], 1.0)

make_body3() = Body(SA[0.0, 0.0, 0.0], SA[-0.93240737, -0.86473146, 0.0], 1.0)

make_figure8_system() = System([make_body1(), make_body2(), make_body3()])

system = make_figure8_system()

pos = simulate(system, Leapfrog(), 0.001, 10000)


function animate_trajectories(positions)
    body1 = Observable(Point2f0(positions[1][1][1], positions[1][1][2]))
    body1_hist = Observable(Point2f0[(positions[1][1][1], positions[1][1][2])])
    body2 = Observable(Point2f0(positions[2][1][1], positions[2][1][2]))
    body2_hist = Observable(Point2f0[(positions[2][1][1], positions[2][1][2])])
    body3 = Observable(Point2f0(positions[3][1][1], positions[3][1][2]))
    body3_hist = Observable(Point2f0[(positions[3][1][1], positions[3][1][2])])

    fig = Figure()
    ax = Axis(fig[1, 1], title="Figure-8 Three-Body Orbit")
    hidespines!(ax, :t, :r)

    ax.xgridvisible = false
    ax.ygridvisible = false

    scatter!(ax, body1; color=:red)
    lines!(ax, body1_hist; color=:red)
    scatter!(ax, body2; color=:green)
    lines!(ax, body2_hist; color=:green)
    scatter!(ax, body3; color=:purple)
    lines!(ax, body3_hist; color=:purple)

    xs = Float64[]
    ys = Float64[]
    for traj ∈ positions, p ∈ traj
        push!(xs, p[1])
        push!(ys, p[2])
    end
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    limits!(ax, xmin - 0.05, xmax + 0.05, ymin - 0.05, ymax + 0.05)

    steps = 1:40:6326

    record(fig, "./docs/images/figure8.mp4", steps;
        framerate=30) do step
        new_point1 = Point2f0(positions[1][step][1], positions[1][step][2])
        new_point2 = Point2f0(positions[2][step][1], positions[2][step][2])
        new_point3 = Point2f0(positions[3][step][1], positions[3][step][2])
        body1[] = new_point1
        body2[] = new_point2
        body3[] = new_point3

        body1_hist[] = push!(body1_hist[], new_point1)
        body2_hist[] = push!(body2_hist[], new_point2)
        body3_hist[] = push!(body3_hist[], new_point3)
    end
end

animate_trajectories(pos)