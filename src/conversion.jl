# Dispersion Gradient
struct DispersionGradient <: Equivalence
    dispersion::Dispersion
    gradient::Union{Function,Nothing}
    gradient_inverse::Union{Function,Nothing}

    function DispersionGradient(dispersion::Dispersion;
        gradient::Union{Function,Nothing}=nothing,
        gradient_inverse::Union{Function,Nothing}=nothing
    )
        return new(dispersion, gradient, gradient_inverse)
    end
end

function DispersionGradient(;
    dispersion::Union{Function,Nothing}=nothing,
    dispersion_inverse::Union{Function,Nothing}=nothing,
    gradient::Union{Function,Nothing}=nothing,
    gradient_inverse::Union{Function,Nothing}=nothing
)
    dispersion = Dispersion(; dispersion, dispersion_inverse)
    return DispersionGradient(dispersion; gradient, gradient_inverse)
end

_frequency_types() = [
    Time, Frequency, AngularPeriod, AngularVelocity,
    Length, Wavenumber, AngularWavelength, AngularWavenumber
]

for T1 ∈ _frequency_types(), T2 ∈ _frequency_types()
    @eval begin
        function edconvert(d::dimtype($T1), x::$T2, dispersion::DispersionGradient)
            return edconvert(d, x, dispersion.dispersion)
        end
    end
end

# uconvert Spectrum
function uconvert(
    uq::Units,
    u1::Units,
    u2::Units,
    x::Spectrum,
    dispersion::DispersionGradient=DispersionGradient();
)
    # axis 1
    if isdirection(x.axis1)
        axis1 = uconvert.(u1, x.axis1)
        g1 = ones(length(x.axis1))
    else
        axis1 = uconvert.(u1, x.axis1, dispersion)
        g1 = _get_gradients(u1, x.axis1, dispersion)
    end
    # axis 2
    if isdirection(x.axis2)
        axis2 = uconvert.(u2, x.axis2)
        g2 = ones(length(x.axis2))
    else
        axis2 = uconvert.(u2, x.axis2, dispersion)
        g2 = _get_gradients(u2, x.axis2, dispersion)
    end
    # data
    data = uconvert.(uq / (u1 * u2), x.data ./ g1 ./ (g2'))
    return Spectrum(data, axis1, axis2)
end

uconvert(uq::Units, x::Spectrum) = uconvert(uq, unit(x, :axis1), unit(x, :axis2), x)

function uconvert(u::Units, s::Symbol, x::Spectrum, dispersion::DispersionGradient=DispersionGradient())
    uq, u1, u2 = unit(x, :integral), unit(x, :axis1), unit(x, :axis2)
    if s==x.axesnames[1] || s==:axis1
        u1 = u
    elseif s==x.axesnames[2] || s==:axis2
        u2 = u
    elseif s==:integral
        uq = u
    elseif s==:spectrum
        throw(ArgumentError(
            "Provide the units of the integral quantity and axes separately."
        ))
    else
        throw(ArgumentError("Invalid quantity."))
    end
    return uconvert(uq, u1, u2, x, dispersion)
end


# uconvert OmnidirectionalSpectrum
function uconvert(
    uq::Units,
    uax::Units,
    x::OmnidirectionalSpectrum,
    dispersion::DispersionGradient=DispersionGradient();
)
    # axis
    if isdirection(x.axis)
        axis = uconvert.(u1, x.axis)
        g1 = ones(length(x.axis))
    else
        axis = uconvert.(uax, x.axis, dispersion)
        g1 = _get_gradients(uax, x.axis, dispersion)
    end
    # data
    data = uconvert.(uq / uax, x.data ./ g1)
    return OmnidirectionalSpectrum(data, axis)
end

uconvert(uq::Units, x::OmnidirectionalSpectrum) = uconvert(uq, unit(x, :axis), x)

function uconvert(
    u::Units,
    s::Symbol,
    x::OmnidirectionalSpectrum,
    dispersion::DispersionGradient=DispersionGradient()
)
    uq, uax = unit(x, :integral), unit(x, :axis)
    if s == x.axisname || s == :axis
        uax = u
    elseif s == :integral
        uq = u
    elseif s == :spectrum
        throw(ArgumentError(
            "Provide the units of the integral quantity and axis separately."
        ))
    else
        throw(ArgumentError("Invalid quantity."))
    end
    return uconvert(uq, uax, x, dispersion)
end


# support functions
_derivative() = Dict(
    # 0 - temporal
    (𝐓, 𝐓) => (x -> 1),
    (𝐓^-1, 𝐓^-1) => (x -> 1),
    (𝐀^-1 * 𝐓, 𝐀^-1 * 𝐓) => (x -> 1),
    (𝐀 * 𝐓^-1, 𝐀 * 𝐓^-1) => (x -> 1),
    # 0 - spatial
    (𝐋, 𝐋) => (x -> 1),
    (𝐋^-1, 𝐋^-1) => (x -> 1),
    (𝐀^-1 * 𝐋, 𝐀^-1 * 𝐋) => (x -> 1),
    (𝐀 * 𝐋^-1, 𝐀 * 𝐋^-1) => (x -> 1),
    # 1 - temporal
    (𝐓^-1, 𝐓) => (x -> -1 / x^2),
    (𝐓, 𝐓^-1) => (x -> -1 / x^2),
    (𝐓^-1, 𝐀 * 𝐓^-1) => (x -> 2π * rad),
    (𝐀 * 𝐓^-1, 𝐓^-1) => (x -> 1 / (2π * rad)),
    (𝐓, 𝐓 * 𝐀^-1) => (x -> 1 / (2π * rad)),
    (𝐓 * 𝐀^-1, 𝐓) => (x -> 2π * rad),
    (𝐀 * 𝐓^-1, 𝐓 * 𝐀^-1) => (x -> -1 / x^2),
    (𝐓 * 𝐀^-1, 𝐀 * 𝐓^-1) => (x -> -1 / x^2),
    # 1 - spatial
    (𝐋^-1, 𝐋) => (x -> -1 / x^2),
    (𝐋, 𝐋^-1) => (x -> -1 / x^2),
    (𝐋^-1, 𝐀 * 𝐋^-1) => (x -> 2π * rad),
    (𝐀 * 𝐋^-1, 𝐋^-1) => (x -> 1 / (2π * rad)),
    (𝐋, 𝐋 * 𝐀^-1) => (x -> 1 / (2π * rad)),
    (𝐋 * 𝐀^-1, 𝐋) => (x -> 2π * rad),
    (𝐀 * 𝐋^-1, 𝐋 * 𝐀^-1) => (x -> -1 / x^2),
    (𝐋 * 𝐀^-1, 𝐀 * 𝐋^-1) => (x -> -1 / x^2),
    # 2 - temporal
    (𝐓^-1, 𝐓 * 𝐀^-1) => (x -> -1 / (2π * rad * x^2)),
    (𝐓 * 𝐀^-1, 𝐓^-1) => (x -> -1 / (2π * rad * x^2)),
    (𝐓^-1 * 𝐀, 𝐓) => (x -> -2π * rad / x^2),
    (𝐓, 𝐓^-1 * 𝐀) => (x -> -2π * rad / x^2),
    # 2 - spatial
    (𝐋^-1, 𝐋 * 𝐀^-1) => (x -> -1 / (2π * rad * x^2)),
    (𝐋 * 𝐀^-1, 𝐋^-1) => (x -> -1 / (2π * rad * x^2)),
    (𝐋^-1 * 𝐀, 𝐋) => (x -> -2π * rad / x^2),
    (𝐋, 𝐋^-1 * 𝐀) => (x -> -2π * rad / x^2),
)

function _get_gradients(u::Units, axis::AbstractVector, dispersion::DispersionGradient)
    from_dim, to_dim = axesinfo(axis)[2], dimension(u)
    if istemporal(from_dim) && isspatial(to_dim)
        ax_t = uconvert.(rad/s, axis, dispersion)
        ax_s = uconvert.(rad / m, axis, dispersion)
        grad_1 = _derivative()[from_dim, 𝐀/𝐓]
        grad_2 = ax -> dispersion.gradient_inverse.(ax)
        grad_3 = _derivative()[𝐀/𝐋, to_dim]
        grad = grad_1.(axis) .* grad_2.(ax_t) .* grad_3.(ax_s)
    elseif isspatial(from_dim) && istemporal(to_dim)
        nothing
        ax_s = uconvert.(rad / m, axis, dispersion)
        ax_t = uconvert.(rad / s, axis, dispersion)
        grad_1 = _derivative()[from_dim, 𝐀/𝐋]
        grad_2 = ax -> dispersion.gradient.(ax)
        grad_3 = _derivative()[𝐀/𝐓, to_dim]
        grad = grad_1.(axis) .* grad_2.(ax_s) .* grad_3.(ax_t)
    else
        grad = (_derivative()[from_dim, to_dim]).(axis)
    end
    return uconvert.(u/unit(eltype(axis)), abs.(grad))
end
