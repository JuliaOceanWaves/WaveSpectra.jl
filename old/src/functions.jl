
# axisnames
axisnames(S::Spectrum) = (axistype(S.axis1), axistype(S.axis2))
axisnames(S::OmniSpectrum) = axistype(S.axis)

# units
function unit(s::Spectrum, quantity::Symbol) :: Units
    return (
        quantity==:axis1 ? unit(eltype(s.axis1)) :
        quantity==:axis2 ? unit(eltype(s.axis2)) :
        quantity==axistype(s.axis1) ? unit(eltype(s.axis1)) :
        quantity==axistype(s.axis2) ? unit(eltype(s.axis2)) :
        quantity==:integral ? unit(eltype(s))*unit(eltype(s.axis1))*unit(eltype(s.axis2)) :
        quantity==:spectrum ? unit(eltype(s)) :
        error("Unknown `quantity`."))
end

unit(s::Spectrum) = unit(s, :spectrum)

function unit(s::OmniSpectrum, quantity::Symbol) :: Units
    return (
        quantity==:axis ? unit(eltype(s.axis)) :
        quantity==axistype(s.axis) ? unit(eltype(s.axis)) :
        quantity==:integral ? unit(eltype(s))*unit(eltype(s.axis)) :
        quantity==:spectrum ? unit(eltype(s)) :
        error("Unknown `quantity`."))
end

unit(s::OmniSpectrum) = unit(s, :spectrum)

# integrals
integrate(S::Spectrum, method::IntegrationMethod=Trapezoidal()) = integrate((S.axis1, S.axis2), S, method)

function integrate(S::Spectrum, ax::Symbol, method::IntegrationMethod=Trapezoidal())
    if (ax==:axis1) || (ax==axisnames(S)[1])
        axis = S.axis2
        data = integrate(S.axis1, S; dims=2)
    elseif (ax==:axis2) || (ax==axisnames(S)[2])
        axis = S.axis1
        data = integrate(S.axis2, S; dims=1)
    else
        error("Unknown axis.")
    end
    return OmniSpectrum(data, axis)
end

integrate(S::OmniSpectrum, method::IntegrationMethod=Trapezoidal()) = integrate(S.axis, S, method)

# change of variable and or units
function convertaxis(u::Units, S::OmniSpectrum; dispersion::Union{Nothing, Equivalence})
    type_org = axistype(S.axis)[1][1]
    type_new = axistype([1*u, ])[1][1]
    if type_org==type_new==:direction
        data = convert.(u, S.data)
        axis = convert.(u, S.axis)
    elseif type_org==type_new
        data = convert.(u, S.data, Periodic())
        axis = convert.(u, S.axis, Periodic())
    elseif type_org==:direction || type_new==:direction
        error("Cannot convert between frequency and direction.")
    elseif isnothing(dispersion)
        error("Converting between temporal and spatial quantities requires a dispersion relation.")
    else
        data =
        axis =
    end
    return OmniSpectrum(data, axis)
end

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
