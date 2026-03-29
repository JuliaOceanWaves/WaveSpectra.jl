"""
Dispersion relations for linear interfacial waves.

Currently only gravity wave dispersion relations are implemented.
The intention is to grow this library.

See also: https://en.wikipedia.org/wiki/Dispersion_(water_waves)
"""
module DispersionRelations

using ..WaveSpectra: Dispersion, g, m, rad, s, θ₀
using Unitful: kg
using Roots: ZeroProblem, solve

"""
    gravitywaves_deepwater()

Deep-water gravity-wave linear dispersion model based on the dispersion
relation: ``ω(k) = √(g k)``.
"""
function gravitywaves_deepwater()
    dispersion = (k -> √(g * k * θ₀))
    dispersion_inverse = (ω -> ω^2 / (g * θ₀))
    gradient = (k -> √(g * θ₀ / k) / 2)
    return Dispersion(; dispersion, dispersion_inverse, gradient)
end

"""
    gravitywaves_shallowwater(waterdepth::Length)

Shallow-water gravity-wave linear dispersion model based on the dispersion
relation: ``ω(k) = k √(g h)``.

The `waterdepth` must have units of length.
Assumes zero air density.
"""
function gravitywaves_shallowwater(waterdepth)
    h = waterdepth
    dispersion = (k -> k * √(g * h))
    dispersion_inverse = (ω -> ω / √(g * h))
    gradient = (k -> √(g * h))
    return Dispersion(; dispersion, dispersion_inverse, gradient)
end

"""
    gravitywaves(waterdepth::Length)

Finite-depth linear gravity-wave dispersion model based on the dispersion
relation: ``ω(k) = √(g k tanh(k h))``.

The `waterdepth` must have units of length.
Assumes zero air density.
"""
function gravitywaves(waterdepth)
    h = waterdepth
    dispersion = (k -> √(k * θ₀ * g * tanh(k * h / θ₀)))
    function dispersion_inverse(ω)
        f = (k -> k - ω^2 / (g * θ₀ * tanh(k * h / θ₀)))
        k₀ = (2π) * rad / m  # initial guess: 1m wavelength
        return solve(ZeroProblem(f, k₀))
    end
    function gradient(k)
        kh = k * h / θ₀
        cp = √((g * θ₀ / k) * tanh(kh))
        return 0.5 * cp * (1 + 2kh / sinh(2kh))
    end
    return Dispersion(; dispersion, dispersion_inverse, gradient)
end

end
