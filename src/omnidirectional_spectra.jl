const _frequency_dims = [𝐓, 𝐓^-1, 𝐓 * 𝐀^-1, 𝐀 * 𝐓^-1, 𝐋, 𝐋^-1, 𝐋 * 𝐀^-1, 𝐀 * 𝐋^-1, NoDims]
const _temporal_frequency_types = [Time, Frequency, AngularPeriod, AngularVelocity]
const _spatial_frequency_types = [Length, Wavenumber, AngularWavelength, AngularWavenumber]
const _frequency_types = vcat(_temporal_frequency_types, _spatial_frequency_types)
const _Temporal = Union{Time,Frequency,AngularPeriod,AngularVelocity}
const _Spatial = Union{Length,Wavenumber,AngularWavelength,AngularWavenumber}
const _TF_spatial = 1.0 * rad / m
const _TF_temporal = 1.0 * rad / s

# dispersion gradient struct
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
    da_dispersion = Dispersion(dispersion, dispersion_inverse)
    return DispersionGradient(da_dispersion; gradient, gradient_inverse)
end

for T1 ∈ _frequency_types, T2 ∈ _frequency_types
    @eval begin
        function UnitfulEquivalences.edconvert(d::dimtype($T1), x::$T2, dispersion::DispersionGradient)
            return edconvert(d, x, dispersion.dispersion)
        end
    end
end

# dispersion gradient func
const _grad_1 = Dict(
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
    (𝐓^-1, 𝐓) => (x -> -1 / x^2), # TODO Check
    (𝐓, 𝐓^-1) => (x -> -1 / x^2), # TODO Check
    (𝐓^-1, 𝐀 * 𝐓^-1) => (x -> 2π * rad),
    (𝐀 * 𝐓^-1, 𝐓^-1) => (x -> 1 / (2π * rad)),
    (𝐓, 𝐓 * 𝐀^-1) => (x -> 1 / (2π * rad)),
    (𝐓 * 𝐀^-1, 𝐓) => (x -> 2π * rad),
    (𝐀 * 𝐓^-1, 𝐓 * 𝐀^-1) => (x -> -1 / x^2), # TODO Check
    (𝐓 * 𝐀^-1, 𝐀 * 𝐓^-1) => (x -> -1 / x^2), # TODO Check
    # 1 - spatial
    (𝐋^-1, 𝐋) => (x -> -1 / x^2), # TODO Check
    (𝐋, 𝐋^-1) => (x -> -1 / x^2), # TODO Check
    (𝐋^-1, 𝐀 * 𝐋^-1) => (x -> 2π * rad),
    (𝐀 * 𝐋^-1, 𝐋^-1) => (x -> 1 / (2π * rad)),
    (𝐋, 𝐋 * 𝐀^-1) => (x -> 1 / (2π * rad)),
    (𝐋 * 𝐀^-1, 𝐋) => (x -> 2π * rad),
    (𝐀 * 𝐋^-1, 𝐋 * 𝐀^-1) => (x -> -1 / x^2), # TODO Check
    (𝐋 * 𝐀^-1, 𝐀 * 𝐋^-1) => (x -> -1 / x^2), # TODO Check
)

_get_grad(dims_from::Dimensions, dims_to::Dimensions) = _grad_1[dims_from, dims_to]
_firstState() = 1

include("omnidirectional_spectra/discrete.jl")
include("omnidirectional_spectra/continuous.jl")
include("omnidirectional_spectra/wave_time_series.jl")
include("omnidirectional_spectra/ocean_waves.jl")


"""
    frequency_unit(::[Discrete]OmnidirectionalSpectrum{TS, TF})

Return the dimension of the expected frequency vector.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1.0u"Hz", 5.0u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v,f);

julia> s2 = OmnidirectionalSpectrum(s1);

julia> unit(s1) == unit(s2) == u"Hz"
true

```
"""
function frequency_unit end

"""
    frequency_dimension(::[Discrete]OmnidirectionalSpectrum{TS, TF})

Return the dimension of the expected frequency vector.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> f=v=range(1.0u"Hz", 5.0u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v,f);

julia> s2 = OmnidirectionalSpectrum(s1);

julia> dimension(s1) == dimension(s2) == dimension(u"Hz")
true

```
"""
function frequency_dimension end

"""
    quantity(::[Discrete]OmnidirectionalSpectrum{TS, TF})

Return the dimensions and units of the product between spectra and frequency.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> f=v=range(1.0u"Hz", 5.0u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v,f);

julia> s2 = OmnidirectionalSpectrum(s1);

julia> quantity(s1)
(𝐓⁻², Hz²)

julia> quantity(s2)
(𝐓⁻², Hz²)

```
"""
function quantity end

function DiscreteOmnidirectionalSpectrum(timeseries::WaveTimeSeries{TS}, frequency::AbstractRange{<:Quantity}) where {TS}
    n_samples = size(timeseries.value, 1)÷2+1
    ϕ = im .*randn(n_samples)
    rfft_temp = (rfft(ustrip(timeseries.value)))
    Δfreq = step(frequency)
    value = @. abs((1/n_samples^2)*(rfft_temp.*unit(TS) * exp(ϕ))^2 / (2*Δfreq))
    return DiscreteOmnidirectionalSpectrum(value, frequency)
end

function DiscreteOmnidirectionalSpectrum(timeseries::WaveTimeSeries{TS}, frequency::AbstractVector{<:Quantity}) where {TS}
    n_samples = size(timeseries.value, 1)÷2+1
    ϕ = im .*randn(n_samples)
    rfft_temp = (rfft(ustrip(timeseries.value)))
    Δfreq = isapprox(0.0, frequency[begin]) ? frequency[2] : frequency[begin]
    value = @. abs((1/n_samples^2)*(rfft_temp.*unit(TS) * exp(ϕ))^2 / (2*Δfreq))
    return DiscreteOmnidirectionalSpectrum(value, frequency)
end