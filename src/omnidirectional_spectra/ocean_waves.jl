
# Parametric spectra
function pierson_moskowitz_spectrum(significant_waveheight::Length, energy_period::Time)
    fₚ = uconvert(Hz, 0.858/energy_period)
    hₛ = significant_waveheight

    function spectrum(f::Frequency)
        if f==0Hz
            value = 0m^2/Hz
        else
            value = (hₛ^2 / 4) * (1.057 * fₚ)^4 * (f)^(-5) * exp((-5 / 4) * (fₚ / f)^4)
        end
        return value
    end
    return OmnidirectionalSpectrum(spectrum)
end

# Deep water dispersion relation
deepwater = Dispersion(
    dispersion=(k -> √(g * k * θ₀)),
    dispersion_inverse=(ω -> ω^2 / (g * θ₀))
);

const deepwater_gradient = DispersionGradient(
    deepwater;
    gradient=(k -> √(gθ₀ / 4k)),
    gradient_inverse=(ω -> 2ω / (g * θ₀))
);

# Spread function


# Statistics (some useful ones, not comprehensive)
"""
    energy_period(
        spectrum::OmnidirectionalSpectrum{TS,TF},
        f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...
    ) where {TS,TF<:Frequency}

Calculate the energy period of an omnidirectional spectra between a range if given.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"));

julia> energy_period(s, 1u"Hz", 5u"Hz")
0.3333333333333333 s

```
"""
function energy_period(
        spectrum::OmnidirectionalSpectrum{TS,TF},
        f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...
    ) where {TS,TF<:Frequency}
    m₋₁ = spectral_moment(spectrum, -1, f_begin, f_end; alg, kwargs...)
    m₀ = spectral_moment(spectrum, 0, f_begin, f_end; alg, kwargs...)
    return m₋₁ / m₀
end

"""
    energy_period(spectrum::OmnidirectionalSpectrum{TS,TF},
        dispersion::Dispersion=Dispersion(), f_begin::Union{Quantity,Nothing}=nothing,
        f_end::Union{Quantity,Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL(),
        kwargs...
    ) where {TS,TF}

Calculate the energy period of an omnidirectional spectra in Hz.

# Example
```jldoctest
julia> using WaveSpectra, DimensionfulAngles

julia> using Unitful:@u_str

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz^-1"));

julia> energy_period(s, Dispersion(), 1u"Hz", 5u"Hz")
0.688888888888889 s

```
"""
function energy_period(spectrum::OmnidirectionalSpectrum{TS,TF},
    dispersion::Dispersion=Dispersion(), f_begin::Union{Quantity,Nothing}=nothing,
    f_end::Union{Quantity,Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL(),
    kwargs...) where {TS,TF}
spectrum_hz = convert_frequency(spectrum, typeof(one(TF) * Hz), dispersion)
return energy_period(spectrum_hz, f_begin, f_end; alg, kwargs...)
end

"""
    energy_period(spectrum::DiscreteOmnidirectionalSpectrum; alg::AbstractIntegralAlgorithm=TrapezoidalRule())

Calculate the energy period of a discrete spectra.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> energy_period(s)
0.3333333333333333 s

```
"""
function energy_period(spectrum::DiscreteOmnidirectionalSpectrum;
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(spectrum, -1; alg)
    m_0 = spectral_moment(spectrum, 0; alg)
    return m_n1 ./ m_0
end

"""
    significant_waveheight(spectrum::OmnidirectionalSpectrum{TS,TF},
        f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...
    ) where {TS,TF}

Calculate the significant waveheight of an omnidirectional spectra.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> significant_waveheight(s)
13.856406460551018 s⁻¹

```
"""
function significant_waveheight(spectrum::OmnidirectionalSpectrum{TS,TF},
    f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
    alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...) where {TS,TF}
    @assert quantity(spectrum)[1] == 𝐋^2
    m₀ = spectral_moment(spectrum, 0, f_begin, f_end; alg, kwargs...)
    return 4√m₀
end

"""
    significant_waveheight(spectrum::DiscreteOmnidirectionalSpectrum; alg::AbstractIntegralAlgorithm=TrapezoidalRule())

Calculate the significant waveheight of a discrete spectra.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> significant_waveheight(s)
13.856406460551018 s⁻¹

```
"""
function significant_waveheight(spectrum::DiscreteOmnidirectionalSpectrum;
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_0 = spectral_moment(spectrum, 0; alg)
    return @. 4*√m_0
end

function steepness(spectrum::OmnidirectionalSpectrum{TS,TF}, dispersion::Dispersion,
    f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
    alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...) where {TS,TF}
    hₛ = significant_waveheight(spectrum, args...; alg, kwargs...)
    tₑ = energy_period(spectrum, f_begin, f_end; alg, kwargs...)
    λₑ = uconvert(m, tₑ, dispersion)
    return hₛ / λₑ
end

# normalize/shape
function normalize(spectrum::OmnidirectionalSpectrum{TS,TF},
        f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(), dispersion::Dispersion=Dispersion(),
        kwargs...) where {TS, TF}
    spectrum_hz = convert_frequency(spectrum, typeof(one(TF) * Hz), dispersion)
    hₛ = significant_waveheight(spectrum_hz, f_begin, f_end; alg, kwargs...)
    tₑ = energy_period(spectrum_hz, f_begin, f_end; alg, kwargs...)
    func(f::DimensionlessQuantity) = uconvert(∅, spectrum(f / tₑ) / (hₛ^2 * tₑ))
    return OmnidirectionalSpectrum(func, typeof(1.0∅))
end

function scale(spectrum::OmnidirectionalSpectrum{TS,TF}, significant_waveheight::Length,
        energy_period::Time) where {TS<:DimensionlessQuantity, TF<:DimensionlessQuantity}
    hₛ, tₑ = (significant_waveheight, energy_period)
    func(f::Frequency) = uconvert(m^2/Hz, spectrum.func(f*tₑ) * (hₛ^2 * tₑ))
    TFnew = typeof(one(typeof(tₑ))*Hz)
    return OmnidirectionalSpectrum(func, TFnew)
end

# Slope spectrum #TODO: check
function slope_spectrum(spectrum::OmnidirectionalSpectrum{TS,TF}) where {TS,TF<:AngularWavenumber}
    func(k::AngularWavenumber) = k^2 * spectrum(k)
    return OmnidirectionalSpectrum(func, TF)
end

function slope_spectrum(spectrum::OmnidirectionalSpectrum{TS,TF},
    dispersion::Dispersion=Dispersion()) where {TS,TF}
    spectrum_rad_m = convert_frequency(spectrum, typeof(one(TF) * rad / m), dispersion)
    slope_rad_m = slope_spectrum(spectrum_rad_m)
    return convert_frequency(slope_rad_m, TF, dispersion)
end

function phase_velocity(frequency::Quantity, dispersion::Dispersion)
    wavelength = uconvert(m, frequency, dispersion)
    period = uconvert(s, frequency, dispersion)
    return wavelength / period
end

function phase_velocity(frequency::Quantity, dispersion::DispersionGradient)
    return phase_velocity(frequency, dispersion.dispersion)
end
