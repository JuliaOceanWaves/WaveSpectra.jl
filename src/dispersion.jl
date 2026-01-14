
deepwater = DispersionGradient(;
    dispersion = (k -> √(g * k * θ₀)),
    dispersion_inverse = (ω -> ω^2 / (g * θ₀)),
    gradient = (k -> √(g/k)/2),
    gradient_inverse = (ω -> 2ω/g),
)
