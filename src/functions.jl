
# units
function unit(x::Spectrum, quantity::Symbol)::Units
    ux, u1, u2 = unit(eltype(x)), unit(eltype(x.axis1)), unit(eltype(x.axis2))
    (quantity == :axis1) && return u1
    (quantity == :axis2) && return u2
    (quantity == x.axesnames[1]) && return u1
    (quantity == x.axesnames[2]) && return u2
    (quantity == :integral) && return ux * u1 * u2
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::Spectrum) = unit(x, :spectrum)

function unit(x::OmnidirectionalSpectrum, quantity::Symbol)::Units
    ux, ua = unit(eltype(x)), unit(eltype(x.axis))
    (quantity == :axis) && return ua
    (quantity == axestypes(x.axis)) && return ua
    (quantity == :integral) && return ux * ua
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::OmnidirectionalSpectrum) = unit(x, :spectrum)


# utilities
AxisArrays.axes(x::Spectrum) = (x.axis1, x.axis2)
AxisArrays.axes(x::OmnidirectionalSpectrum) = (x.axis,)

coordinates(x::Spectrum) = x.coordinates
axesnames(x::Spectrum) = x.axesnames
axesnames(x::OmnidirectionalSpectrum) = x.axisname

ispolar(x::Spectrum) = (x.coordinates == :polar)
iscartesian(x::Spectrum) = (x.coordinates == :cartesian)




# # change of variable and or units
# function convertaxis(u::Units, S::OmniSpectrum; dispersion::Union{Nothing, Equivalence})
#     type_org = axistypes(S.axis)[1][1]
#     type_new = axistypes([1*u, ])[1][1]
#     if type_org==type_new==:direction
#         data = convert.(u, S.data)
#         axis = convert.(u, S.axis)
#     elseif type_org==type_new
#         data = convert.(u, S.data, Periodic())
#         axis = convert.(u, S.axis, Periodic())
#     elseif type_org==:direction || type_new==:direction
#         error("Cannot convert between frequency and direction.")
#     elseif isnothing(dispersion)
#         error("Converting between temporal and spatial quantities requires a dispersion relation.")
#     # else
#     #     data =
#     #     axis =
#     # end
#     return OmniSpectrum(data, axis)
# end

# derivative(f, x::Quantity) = derivative(f, ustrip(x)) * unit(f(x)) / unit(x)

# function derivative(f, X::AbstractArray{<:Quantity})
#     _check_typeconsistency(X)
#     return derivative(f, ustrip.(X)) * unit(f(X[0])) / unit(X[0])
# end

# changeaxis(f, S::OmniSpectrum) = OmniSpectrum(S.data / derivative(f, S.axis), f(S.axis))



# S(ω, θ) <=> S(ω), D(ω, θ) # ω, f, T, k, λ, ...


# Quality check: integral D(ω, θ) ≈ 1, θ∈[0°, ..., 360°)


# functions omnispectra: moments, Hs, Tp, etc etc from MHKiT


# time-series <=> Spectrum
# wave elevation type?
# non-even spacing?
# elevation variance spectrum + phase information + randomization of spectra
