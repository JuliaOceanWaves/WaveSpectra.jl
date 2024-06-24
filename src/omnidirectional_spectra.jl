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

_get_grad(dims_from::Dimensions, dims_to::Dimensions) = _grad_1[dims_from, dims_to]

include("omnidirectional_spectra/discrete.jl")
include("omnidirectional_spectra/continuous.jl")
include("omnidirectional_spectra/ocean_waves.jl")
