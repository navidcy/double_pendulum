using DynamicalSystems
using Plots

ENV["GKSwstype"] = "100"

include("double_pendulum_eom.jl")

noforcing(t) = 0

function make_a_movie(θ₁₀ = 0.25, ω₁₀ = 0.2, θ₂₀ = 0.0, ω₂₀ = 0.4;
                      tfinal = 10, dt = 0.01, Δt = 10dt,
                      filename = "movie.mp4",
                      g = 9.81, ℓ₁ = 1.0, ℓ₂ = 1.0, m₁ = 1.0, m₂ = 1.0,
                      γ = 0.0, forcing=noforcing,
                      anglelimit = 30,
                      dpi = 200, fps = 36
                      )

    t = 0:dt:tfinal
    
    params = (g, ℓ₁, ℓ₂, m₁, m₂, γ, forcing)
    u₀ = [θ₁₀, ω₁₀, θ₂₀, ω₂₀]
    
    ds = ContinuousDynamicalSystem(double_pendulum, u₀, params)

    solution = trajectory(ds, tfinal, u₀, dt=dt)

    θ₁_full = solution[:, 1]
    θ₂_full = solution[:, 3]

    global θ₁, ω₁, θ₂, ω₂ = θ₁₀, ω₁₀, θ₂₀, ω₂₀
    
    anim = @animate for j in 0:Int(tfinal/Δt)
        local u = [θ₁, ω₁, θ₂, ω₂]
        
        time = j * Δt
        
        local ds = ContinuousDynamicalSystem(double_pendulum, u, params; t0=time)
        
        local tr = j==0 ? trajectory(ds, 0, u) : trajectory(ds, Δt, u, dt=dt)
        
        global θ₁ = tr[end, 1]
        global ω₁ = tr[end, 2]
        global θ₂ = tr[end, 3]
        global ω₂ = tr[end, 4]
        
        # convert to Cartesian coords to plot
        local x₁ =  ℓ₁ * sin(θ₁)
        local y₁ = -ℓ₁ * cos(θ₁)
        local x₂ = x₁ + ℓ₂ * sin(θ₂)
        local y₂ = y₁ - ℓ₂ * cos(θ₂)
        
        local p1 = plot([0, x₁, x₂], [0, y₁, y₂],
                              color = :black,
                          linewidth = 1.5,
                             legend = :none,
                           showaxis = :false,
                               grid = :off,
                             xticks = :none,
                             yticks = :none,
                        aspectratio = 1,
                              xlims = (-2.2, 2.2),
                              ylims = (-2.2, 1.4))
            
        plot!(p1, [0], [0],
              shape = :square, markersize = 5, color = :black, legend = :none)
        plot!(p1, [x₁], [y₁],
              shape = :circle, markersize = 10, color = :red, legend = :none)
        plot!(p1, [x₂], [y₂],
              shape = :hex, markersize = 10, color = :blue, legend = :none)
        
        plot(p1, size=(500, 500), dpi=dpi)
    end

    gif(anim, filename, fps=fps)
end

# Now make movies

# First, no damping, no forcing
tfinal = 50
Δt = 0.025   # time between snapshots in animations
dt = 0.005   # actual dt for time-integration

θ₁₀ = 5π/6
ω₁₀ = 0.2
θ₂₀ = 1.0
ω₂₀ = 1.4

ℓ₁ = 1.1
ℓ₂ = 0.9

filename = "logo.gif"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀;
             ℓ₁ = ℓ₁, ℓ₂ = ℓ₂,
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename, fps=36, dpi=150)
