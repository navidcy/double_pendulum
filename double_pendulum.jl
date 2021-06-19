using DynamicalSystems
using Plots

ENV["GKSwstype"] = "100"

"""
    double_pendulum(u, params, t)
    
```math
\\begin{aligned}
    θ̇₁ &= ω₁ \\\\
    ω̇₁ &= [m₂ ℓ₁ ω₁² \\sin Δθ \\cos Δθ + m₂ g \\sin θ₂ \\cos Δθ +
       m₂ ℓ₂ ω₂² \\sin Δθ - (m₁ + m₂) g \\sin θ₁] / (ℓ₁ Δ) - γ ω₁ + forcing \\\\
    θ̇₂ &= ω₂ \\\\
    ω̇₂ &= [-m₂ ℓ₂ ω₂² \\sin Δθ \\cos Δθ + (m₁ + m₂) g \\sin θ₁ \\cos Δθ -
             (m₁ + m₂) ℓ₁ ω₁² \\sin Δθ - (m₁ + m₂) g \\sin Θ₂] / (ℓ₂ Δ) - γ ω₂
\\end{aligned}
```
where ``Δθ = θ₂ - θ₁`` and ``D = (m₁ + m₂) - m₂ \\cos² Δθ``.
"""
function double_pendulum(u, params, t)
    # unpack parameters and state vector
    g, ℓ₁, ℓ₂, m₁, m₂, γ, forcing = params
    θ₁, ω₁, θ₂, ω₂ = u
    
    # some definitions
    Δθ = θ₂ - θ₁
    D = (m₁ + m₂) - m₂ * cos(Δθ)^2
    
    # equations of motion
    dθ₁ = ω₁
    
    dω₁ = (m₂ * ℓ₁ * ω₁^2 * sin(Δθ) * cos(Δθ) +
               m₂ * g * sin(θ₂) * cos(Δθ) +
               m₂ * ℓ₂ * ω₂^2 * sin(Δθ) -
               (m₁ + m₂) * g * sin(θ₁)) / (ℓ₁ * D) - γ * ω₁ + forcing(t)
    
    dθ₂ = ω₂
    
    dω₂ = (-m₂ * ℓ₂ * ω₂^2 * sin(Δθ) * cos(Δθ) +
               (m₁ + m₂) * g * sin(θ₁) * cos(Δθ) -
               (m₁ + m₂) * ℓ₁ * ω₁^2 * sin(Δθ) -
               (m₁ + m₂) * g * sin(θ₂)) / (ℓ₂ * D) - γ * ω₂
    
    return SVector{4}(dθ₁, dω₁, dθ₂, dω₂)
end

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
    u₀ = (θ₁₀, ω₁₀, θ₂₀, ω₂₀)
    
    ds = ContinuousDynamicalSystem(double_pendulum, u₀, params)

    solution = trajectory(ds, tfinal, (θ₁₀, ω₁₀, θ₂₀, ω₂₀), dt=dt)

    θ₁_full = solution[:, 1]
    θ₂_full = solution[:, 3]

    global θ₁, ω₁, θ₂, ω₂ = θ₁₀, ω₁₀, θ₂₀, ω₂₀
    
    anim = @animate for j in 0:Int(tfinal/Δt)
        local u = (θ₁, ω₁, θ₂, ω₂)
        
        time = j*Δt
        
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
                              xlims = (-1, 1),
                              ylims = (-2.1, 0.1))
            
        plot!(p1, [-1, 1], [0, 0],
              color = :black, linewidth = 2, legend = :none)
        plot!(p1, [0], [0],
              shape = :square, markersize = 5, color = :black, legend = :none)
        plot!(p1, [x₁], [y₁],
              shape = :circle, markersize = 10, color = :red, legend = :none)
        plot!(p1, [x₂], [y₂],
              shape = :hex, markersize = 10, color = :blue, legend = :none)
        
        local p2 = plot(t[1:Int(Δt/dt)*j], 180/π * θ₁_full[1:Int(Δt/dt)*j],
                            color = :red,
                            alpha = 0.5,
                        linewidth = 3,
                           legend = :none,
                           ylabel = "θ₁ [ᵒ]", 
                            ylims = (-anglelimit, anglelimit),
                            xlims = time>10 ? (-0.2, time*1.02) : (-0.2, 10.2))
                  
        local p3 = plot(t[1:Int(Δt/dt)*j], 180/π * θ₂_full[1:Int(Δt/dt)*j],
                            color = :blue,
                            alpha = 0.5,
                        linewidth = 3,
                           legend = :none,
                           xlabel = "time",
                           ylabel = "θ₂ [ᵒ]",
                            ylims = (-anglelimit, anglelimit),
                            xlims = time>10 ? (-0.2, time*1.02) : (-0.2, 10.2))
        
        plot!(p2, [time], [180/π * θ₁], shape=:circle, color=:red)
        plot!(p3, [time], [180/π * θ₂], shape=:hex, color=:blue)


        l = @layout [a{0.45w} grid(2, 1)]

        plot(p1, p2, p3, layout=l, dpi=dpi)
    end

    mp4(anim, filename, fps=fps)
end

# Now make movies

# First, no damping, no forcing
tfinal = 25
Δt = 0.025   # time between snapshots in animations
dt = 0.005   # actual dt for time-integration

θ₁₀ = 0.25
ω₁₀ = 0.2
θ₂₀ = 0 
ω₂₀ = 0.4
filename = "random.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀;
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)

θ₁₀ = 0.25
ω₁₀ = 0.0
θ₂₀ = +sqrt(2) * θ₁₀ 
ω₂₀ = 0
filename = "normalmode1.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀;
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)

θ₁₀ = 0.25
ω₁₀ = 0.0
θ₂₀ = -sqrt(2) * θ₁₀ 
ω₂₀ = 0
filename = "normalmode2.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀;
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)

# Now add damping and forcing and see if we can excite
# the normal modes.

γ=0.4        # damping coefficient

tfinal = 40
Δt = 0.025   # time between snapshots in animations
dt = 0.005   # actual dt for time-integration

θ₁₀ = 0.25
ω₁₀ = 0.2
θ₂₀ = 0 
ω₂₀ = 0.4

ℓ = 1.0
g = 9.81

ξ₁ = sqrt(2 - sqrt(2)) * sqrt(g/ℓ)
ξ₂ = sqrt(2 + sqrt(2)) * sqrt(g/ℓ)
ξ₃ = 1.43sqrt(2 + sqrt(2)) * sqrt(g/ℓ)

φ₁, φ₂, φ₃ = 2π * rand(3)

forcing(t) = 0.3cos(ξ₁ * t + φ₁)
filename = "excite-normalmode1.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀; γ=γ, forcing=forcing, g=g,
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)

forcing(t) = 0.5cos(ξ₂ * t + φ₂)
filename = "excite-normalmode2.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀; γ=γ, forcing=forcing, g=g,
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)

forcing(t) = 0.3cos(ξ₁ * t + φ₁) + 0.5cos(ξ₂ * t + φ₂) + 0.5cos(ξ₃ * t + φ₃)
filename = "excite-random.mp4"
make_a_movie(θ₁₀, ω₁₀, θ₂₀, ω₂₀; γ=γ, forcing=forcing, g=g,
             tfinal = tfinal, dt=dt, Δt=Δt,
             filename = filename)
