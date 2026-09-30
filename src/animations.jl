using CairoMakie

function animate(
    positions;
    output="./out/animation.mp4",
    title="",
    framerate=30,
    show_trails=true,
    limits=nothing,
    hide_decor=false
)
    N = positions.size[1]
    bodies = Vector{Observable{Point2f0}}(undef, N)
    bodies_hist = Vector{Observable{Vector{Point2f0}}}(undef, N)

    for i ∈ 1:N
        bodies[i] = Observable(Point2f0(positions[i][1][1], positions[i][1][2]))
        bodies_hist[i] = Observable(Point2f0[(positions[i][1][1], positions[i][1][2])])
    end

    fig = Figure()
    ax = Axis(fig[1, 1], title=title)
    if hide_decor
        hidespines!(ax)
        hidedecorations!(ax)
    else
        hidespines!(ax, :t, :r)
    end

    ax.xgridvisible = false
    ax.ygridvisible = false

    for i ∈ 1:N
        scatter!(ax, bodies[i])
    end

    if show_trails
        for i ∈ 1:N
            lines!(ax, bodies_hist[i])
        end
    end

    if limits === nothing
        xs = Float64[]
        ys = Float64[]

        for traj ∈ positions, p ∈ traj
            push!(xs, p[1])
            push!(ys, p[2])
        end

        xmin, xmax = extrema(xs)
        ymin, ymax = extrema(ys)

        limits!(ax, xmin - 0.05, xmax + 0.05, ymin - 0.05, ymax + 0.05)
    else
        limits!(ax, limits...)
    end

    steps = 1:(positions[1].size[1])

    mkpath(dirname(output))
    record(fig, output, steps;
        framerate) do step
        for i ∈ 1:N
            new_point = Point2f0(positions[i][step][1], positions[i][step][2])
            bodies[i][] = new_point

            bodies_hist[i][] = push!(bodies_hist[i][], new_point)
        end
    end
end