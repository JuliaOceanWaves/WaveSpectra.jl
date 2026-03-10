# Dispersion Gradient
struct DispersionGradient <: Equivalence
    dispersion::Dispersion
    gradient::Union{Function,Nothing}
    gradient_inverse::Union{Function,Nothing}

    function DispersionGradient(dispersion::Dispersion;
        gradient::Union{Function,Nothing}=nothing,
    )
        function gradient_inverse(ω)
            k = dispersion.dispersion_inverse(ω)
            return inv(gradient(k))
        end
        return new(dispersion, gradient, gradient_inverse)
    end
end

function DispersionGradient(;
    dispersion::Union{Function,Nothing}=nothing,
    dispersion_inverse::Union{Function,Nothing}=nothing,
    gradient::Union{Function,Nothing}=nothing,
)
    dispersion = Dispersion(; dispersion, dispersion_inverse)
    return DispersionGradient(dispersion; gradient)
end

_frequency_types() = [
    Time, Frequency, AngularPeriod, AngularVelocity,
    Length, Wavenumber, AngularWavelength, AngularWavenumber
]
