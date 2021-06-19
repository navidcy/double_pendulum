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
