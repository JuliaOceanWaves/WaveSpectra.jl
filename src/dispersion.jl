deepwater = Dispersion(
    dispersion = (k -> √(g * k * θ₀)),
    dispersion_inverse = (ω -> ω^2 / (g * θ₀))
);
