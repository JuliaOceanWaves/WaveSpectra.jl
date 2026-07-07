# Convert the axes of spectra

"""
    periodic

Periodic equivalence `periodic = DimensionfulAngles.Dispersion()`.
"""
const periodic = Dispersion()

# Spectrum
"""
    uconvert(
        uq::Units,
        u1::Units,
        u2::Units,
        x::AbstractSpectrum,
        dispersion::Dispersion = Dispersion()
    )
    uconvert(uq::Units, x::AbstractSpectrum)
    uconvert(
        u::Units,
        s::Symbol,
        x::AbstractSpectrum,
        dispersion::Dispersion = Dispersion()
    )

Extend `Unitful.uconvert` for directional spectra.

`uconvert(uq, u1, u2, x, dispersion)` converts the integral quantity and both axes of `x`,
updating the spectral-density units so that the integrated quantity is preserved.

`uconvert(uq, x)` converts only the integral quantity and keeps both axes in their current
units.

`uconvert(u, s, x, dispersion)` converts one component selected by `s`, where `s` can be
`:axis1`, `:axis2` (or their axis names), or `:integral`.
Passing `:spectrum` is not supported; the integral quantity and axes units must be
specified explicitly.
"""
function uconvert(
        uq::Units,
        u1::Units,
        u2::Units,
        x::AbstractSpectrum,
        dispersion::Dispersion = Dispersion();
)
    # axis 1
    if isdirection(x.axis1)
        axis1 = uconvert.(u1, x.axis1)
        g1 = ones(length(x.axis1))
        p1 = nothing
    else
        axis1, g1, p1 = _convert_axis(u1, x.axis1, dispersion)
    end
    # axis 2
    if isdirection(x.axis2)
        axis2 = uconvert.(u2, x.axis2)
        g2 = ones(length(x.axis2))
        p2 = nothing
    else
        axis2, g2, p2 = _convert_axis(u2, x.axis2, dispersion)
    end
    data = x.data
    p1 === nothing || (data = data[p1, :])
    p2 === nothing || (data = data[:, p2])
    # data
    data = uconvert.(uq / (u1 * u2), data ./ g1 ./ (g2'))
    return rebuild_superposition(x, data, axis1, axis2)
end

function uconvert(
        u::Units, s::Symbol, x::AbstractSpectrum, dispersion::Dispersion = Dispersion())
    uq, u1, u2 = unit(x, :integral), unit(x, :axis1), unit(x, :axis2)
    if s == x.axesnames[1] || s == :axis1
        u1 = u
    elseif s == x.axesnames[2] || s == :axis2
        u2 = u
    elseif s == :integral
        uq = u
    elseif s == :spectrum
        throw(ArgumentError(
            "Provide the units of the integral quantity and axes separately."
        ))
    else
        throw(ArgumentError("Invalid quantity."))
    end
    return uconvert(uq, u1, u2, x, dispersion)
end

uconvert(uq::Units, x::AbstractSpectrum) = uconvert(uq, :integral, x)

# OmnidirectionalSpectrum
"""
    uconvert(
        uq::Units,
        uax::Units,
        x::AbstractOmnidirectionalSpectrum,
        dispersion::Dispersion = Dispersion()
    )
    uconvert(uq::Units, x::AbstractOmnidirectionalSpectrum)
    uconvert(
        u::Units,
        s::Symbol,
        x::AbstractOmnidirectionalSpectrum,
        dispersion::Dispersion = Dispersion()
    )

Extend `Unitful.uconvert` for omnidirectional spectra.

`uconvert(uq, uax, x, dispersion)` converts the integral quantity and axis of `x`, updating
the spectral-density units so that the integrated quantity is preserved.

`uconvert(uq, x)` converts only the integral quantity and keeps the axis in its current
units.

`uconvert(u, s, x, dispersion)` converts one component selected by `s`, where `s` can be
`:axis` (or its axis name), or `:integral`.
Passing `:spectrum` is not supported; the integral quantity and axis units must be
specified explicitly.
"""
function uconvert(
        uq::Units,
        uax::Units,
        x::AbstractOmnidirectionalSpectrum,
        dispersion::Dispersion = Dispersion();
)
    # axis
    if isdirection(x.axis)
        axis = uconvert.(uax, x.axis)
        g1 = ones(length(x.axis))
        p = nothing
    else
        axis, g1, p = _convert_axis(uax, x.axis, dispersion)
    end
    # data
    data = x.data
    p === nothing || (data = data[p])
    data = uconvert.(uq / uax, data ./ g1)
    return rebuild_superposition(x, data, axis)
end

function uconvert(
        u::Units,
        s::Symbol,
        x::AbstractOmnidirectionalSpectrum,
        dispersion::Dispersion = Dispersion()
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

uconvert(uq::Units, x::AbstractOmnidirectionalSpectrum) = uconvert(uq, :integral, x)

# polar < --- > cartesian
"""
    polar_to_cartesian(x::AbstractSpectrum)

Convert a polar spectrum to cartesian point coordinates.

The result is returned as a 3-column `AxisArray` with one row per sample point:
`(axis_x, axis_y, spectrum)`.
"""
function polar_to_cartesian(x::AbstractSpectrum)
    ispolar(x) || throw(ArgumentError("Spectrum must be in polar coordinates."))
    spectral_axis = x.axis1
    data = x.data
    direction_axis = x.axis2
    ns = length(spectral_axis)
    nd = length(direction_axis)
    n = ns * nd

    axis_x = Vector{typeof(spectral_axis[1] * cos(direction_axis[1]))}(undef, n)
    axis_y = similar(axis_x)
    spectrum = Vector{typeof(data[1, 1] * θ₀ / spectral_axis[1])}(undef, n)

    idx = 1
    @inbounds for j in eachindex(direction_axis)
        cθ = cos(direction_axis[j])
        sθ = sin(direction_axis[j])
        for i in eachindex(spectral_axis)
            kval = spectral_axis[i]
            axis_x[idx] = kval * cθ
            axis_y[idx] = kval * sθ
            spectrum[idx] = data[i, j] * θ₀ / kval
            idx += 1
        end
    end

    name = x.axesnames[1]
    components = [Symbol(string(name) * "_x"), Symbol(string(name) * "_y"), :spectrum]
    return _point_cloud_axisarray(axis_x, axis_y, spectrum, components)
end

"""
    cartesian_to_polar(x::AbstractSpectrum; angle_unit::Units = °)

Convert a cartesian spectrum to polar point coordinates.

Both spectrum axes must have the same spectral-variable type.
The result is returned as a 3-column `AxisArray` with one row per sample point:
`(axis, direction, spectrum)`.
"""
function cartesian_to_polar(x::AbstractSpectrum; angle_unit::Units = °)
    iscartesian(x) || throw(ArgumentError("Spectrum must be in cartesian coordinates."))
    (x.axestypes[1] == x.axestypes[2]) || throw(ArgumentError(
        "Cartesian spectrum axes must have the same spectral-variable type."
    ))

    axis1 = x.axis1
    axis2 = x.axis2
    n1 = length(axis1)
    n2 = length(axis2)
    n = n1 * n2

    spectral_axis = Vector{typeof(sqrt(axis1[1]^2 + axis2[1]^2))}(undef, n)
    direction_axis = Vector{typeof(atan(angle_unit, axis2[1], axis1[1]))}(undef, n)
    spectrum = Vector{typeof(x.data[1, 1] * sqrt(axis1[1]^2 + axis2[1]^2) / θ₀)}(undef, n)

    idx = 1
    @inbounds for j in eachindex(axis2)
        ky = axis2[j]
        for i in eachindex(axis1)
            kx = axis1[i]
            kval = sqrt(kx^2 + ky^2)
            spectral_axis[idx] = kval
            direction_axis[idx] = atan(angle_unit, ky, kx)
            spectrum[idx] = x.data[i, j] * kval / θ₀
            idx += 1
        end
    end

    return _point_cloud_axisarray(
        spectral_axis,
        direction_axis,
        spectrum,
        [x.axestypes[1], :direction, :spectrum]
    )
end

# utilities
function _point_cloud_axisarray(a, b, c, components)
    n = length(a)
    data = Matrix{Any}(undef, n, 3)
    @inbounds for i in 1:n
        data[i, 1] = a[i]
        data[i, 2] = b[i]
        data[i, 3] = c[i]
    end
    return AxisArray(data, Axis{:point}(1:n), Axis{:component}(components))
end

function _convert_axis(u::Units, axis::AbstractVector, dispersion::Dispersion)
    axis_converted = uconvert.(u, axis, dispersion)
    gradients = _get_gradients(u, axis, dispersion)
    mixed_signs = (first(axis) < zero(first(axis))) && (last(axis) > zero(last(axis)))
    f2p = (isfrequency(axis) && isperiod(u)) || (isperiod(axis) && isfrequency(u))
    if f2p && mixed_signs
        perm = sortperm(axis_converted)
        return axis_converted[perm], gradients[perm], perm
    end
    return axis_converted, gradients, nothing
end

function _derivative()
    Dict(
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
        (𝐋, 𝐋^-1 * 𝐀) => (x -> -2π * rad / x^2)
    )
end

function _get_gradients(u::Units, axis::AbstractVector, dispersion::Dispersion)
    from_dim, to_dim = axesinfo(axis)[2], dimension(u)
    if istemporal(from_dim) && isspatial(to_dim)
        ax_t = uconvert.(rad / s, axis, dispersion)
        ax_s = uconvert.(rad / m, axis, dispersion)
        grad_1 = _derivative()[from_dim, 𝐀 / 𝐓]
        grad_2 = ax -> dispersion.gradient_inverse.(ax)
        grad_3 = _derivative()[𝐀 / 𝐋, to_dim]
        grad = grad_1.(axis) .* grad_2.(ax_t) .* grad_3.(ax_s)
    elseif isspatial(from_dim) && istemporal(to_dim)
        ax_s = uconvert.(rad / m, axis, dispersion)
        ax_t = uconvert.(rad / s, axis, dispersion)
        grad_1 = _derivative()[from_dim, 𝐀 / 𝐋]
        grad_2 = ax -> dispersion.gradient.(ax)
        grad_3 = _derivative()[𝐀 / 𝐓, to_dim]
        grad = grad_1.(axis) .* grad_2.(ax_s) .* grad_3.(ax_t)
    else
        grad = (_derivative()[from_dim, to_dim]).(axis)
    end
    return uconvert.(u / unit(eltype(axis)), abs.(grad))
end
