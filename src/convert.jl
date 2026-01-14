
function uconvert(
    uq::Units,
    u1::Units,
    u2::Units,
    x::Spectrum,
    dispersion::Dispersion=Dispersion();
)
    # axis 1
    if isdirection(x.axis1)
        axis1 = uconvert.(u1, x.axis1)
        g1 = (x->1)
    else
        axis1 = uconvert.(u1, x.axis1, dispersion)
        g1 = _derivative()[axesinfo(x.axis1)[2], dimension(u1)]
    end
    # axis 2
    if isdirection(x.axis2)
        axis2 = uconvert.(u2, x.axis2)
        g2 = (x -> 1)
    else
        axis2 = uconvert.(u2, x.axis2, dispersion)
        g2 = _derivative()[axesinfo(x.axis2)[2], dimension(u2)]
    end
    # data
    data = uconvert.(uq / (u1 * u2), x.data ./ abs.(g1.(x.axis1)) ./ abs.((g2.(x.axis2))'))
    return Spectrum(data, axis1, axis2)
end

uconvert(uq::Units, x::Spectrum) = uconvert(uq, unit(x, :axis1), unit(x, :axis2), x)

function uconvert(uq::Units, x::Spectrum, e::Equivalence)
    return uconvert(uq, unit(x, :axis1), unit(x, :axis2), x)
end

function uconvert(u::Units, s::Symbol, x::Spectrum, dispersion::Dispersion=Dispersion())
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
)
