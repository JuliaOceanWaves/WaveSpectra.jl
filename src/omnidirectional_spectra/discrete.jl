# Struct
"""
    DiscreteOmnidirectionalSpectrum(value::AbstractVecOrMat{<:Quantity}, frequency::AbstractVector{<:Quantity}; density::Bool=true)

Create discrete omnidirectional spectrum using two Unitful vectors, or matrix and vector.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v1=f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v1,f);

julia> v2 = ones(typeof(1u"Hz"), 3, 3)
3×3 Matrix{Quantity{Int64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}}:
 1 Hz  1 Hz  1 Hz
 1 Hz  1 Hz  1 Hz
 1 Hz  1 Hz  1 Hz

julia> s2 = DiscreteOmnidirectionalSpectrum(v2, f);

```
"""
struct DiscreteOmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity, D, N} <: AbstractArray{TS, N}
    value::AbstractVecOrMat{<:Quantity}
    frequency::AbstractVector{<:Quantity}
    colnames::Union{AbstractVector{<:AbstractString}, Nothing}
    meta::Any
    function DiscreteOmnidirectionalSpectrum(value::AbstractVecOrMat{<:Quantity}, frequency::AbstractVector{<:Quantity}, colnames=nothing, meta=nothing; density::Bool=true)
        # Parameters
        N=ndims(value) # Don't need to check N if AbstractVecOrMat handles 0 dim and 3+ dims
        D=density
        TF=eltype(frequency)
        TS=eltype(value)
        # Check Parameters
        (dimension(TF) ∉ _frequency_dims) && error("invalid frequency dimensions")
        # Check Arguments
        (size(value)[1] ≠ length(frequency)) && error("'frequency' and first dimension of 'value' array must be same size")
        return new{TS,TF,D,N}(value, frequency)
    end
end

# Constructor
"""
    DiscreteOmnidirectionalSpectrum(func::Function, frequency::AbstractVector{<:Quantity}; density::Bool=true)

Create discrete omnidirectional spectrum using a function and frequency vector.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> func(x) = x
func (generic function with 1 method)

julia> s = DiscreteOmnidirectionalSpectrum(func, f);

```
"""
function DiscreteOmnidirectionalSpectrum(func::Function, frequency::AbstractVector{<:Quantity}; density::Bool=true)
    value = func.(frequency)
    return DiscreteOmnidirectionalSpectrum(value, frequency; density)
end

_firstState() = 1
#Interface
Base.size(spectrum::DiscreteOmnidirectionalSpectrum) = size(spectrum.value)
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, i::Int) = spectrum.value[i]
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, I::Vararg{Int, 2}) = spectrum.value[I...]
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, I) = [spectrum.value[i] for i in I]
Base.length(spectrum::DiscreteOmnidirectionalSpectrum) = length(spectrum.value)
Base.iterate(spectrum::DiscreteOmnidirectionalSpectrum) = isempty(spectrum.value) ? nothing : (spectrum.value[_firstState()], _firstState() + 1)
Base.iterate(spectrum::DiscreteOmnidirectionalSpectrum, state::Int) = (state <= length(spectrum.value)) ? (spectrum.value[state], state+1) : nothing
Base.firstindex(::DiscreteOmnidirectionalSpectrum) = _firstState()
Base.lastindex(spectrum::DiscreteOmnidirectionalSpectrum) = length(spectrum.value)
Base.eltype(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = TS

function Base.show(io::IO, spectrum::DiscreteOmnidirectionalSpectrum{TS,TF,D,N}) where {TS,TF,D,N}
    print(io, "\nPARAMS: {$TS,\n\t $TF,\n\t Density=$D, Value Dims=$N}\nFREQUENCY {$(unit(TF))}, VALUES $(size(spectrum))\t{$(unit(TS))}")
end

function Base.show(io::IO, ::MIME"text/plain", spectrum::DiscreteOmnidirectionalSpectrum{TS,TF,D,N}) where {TS,TF,D,N}
    af = ustrip.(spectrum.frequency)
    if N < 2
        av = ustrip.(spectrum.value)
    else
        lrow_str = "\n\t\t\t⋮   ⋱   ⋮\n\t\t\t"
        erow = [round.(r,digits=3) for r in eachrow(ustrip.(spectrum.value))]
        rows = [length(r) > 2 ? join([r[begin], r[end]], " … ") : join(r, ", ") for r in erow]
        av = length(rows) > 2 ? join([rows[begin], rows[end]], lrow_str) : join(rows, ", \n")
    end
    print(io, "PARAMS: {$TS,\n\t $TF,\n\t Density=$D, Value Dims=$N}\nFREQUENCY\t{$(unit(TF))}: \t$(af)\nVALUES $(size(spectrum))\t{$(unit(TS))}: \t$(av)")
end

# Unitful Interface
Unitful.unit(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = unit(TS)
Unitful.dimension(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = dimension(TS)

# Plots recipes
@recipe function f(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, D}, args...) where {TS, TF, D}
    xlabel --> "frequency"
    ylabel --> (D ? "spectral density" : "discrete (integral) spectrum")
    marker := :auto
    (spectrum.frequency, spectrum.value, args...)
end

"""
    spectral_moment(spectrum::DiscreteOmnidirectionalSpectrum, n::Real=0; alg::AbstractIntegralAlgorithm=TrapezoidalRule())

Calculate the n-th spectral moment of a discrete spectra.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> spectral_moment(s, -1)
4.0 s⁻¹

julia> spectral_moment(s, 0)
12.0 s⁻²

julia> spectral_moment(s, 1)
42.0 s⁻³
```
"""
function spectral_moment(spectrum::DiscreteOmnidirectionalSpectrum, n::Real=0;
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    # There are no keyword arguments used to solve SampledIntegralProblems
    # https://docs.sciml.ai/Integrals/stable/basics/SampledIntegralProblem/
    sol = solve(SampledIntegralProblem(spectrum.value.*(spectrum.frequency .^float(n)), spectrum.frequency; dim=1), alg)
    (sol.retcode != ReturnCode.Success) && error("Solution unsuccessful with code $(sol.retcode)")
    return upreferred.(sol.u)
end

for T1 in [_Temporal, _Spatial]
    @eval begin
        function convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS,TF}, TF_new::T,
                dispersion::Equivalence=deepwater_gradient) where {TS,TF<:$T1, T<:$T1}
            # println("Both Temporal/Spatial")
            grad = _get_grad(dimension(TF), dimension(T))
            new_value = @. upreferred(spectrum.value / grad(spectrum.frequency))
            new_freq = uconvert.(unit(T), spectrum.frequency, dispersion)
            return DiscreteOmnidirectionalSpectrum(new_value, new_freq)
        end
    end
end

function convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS,TF}, TF_new::T,
            dispersion::Equivalence=deepwater_gradient) where {TS,TF<:_Spatial,T<:_Temporal}

    spectrum_int_spatial = convert_frequency(spectrum, _TF_spatial, dispersion)

    grad = dispersion.gradient
    inter_value = @. upreferred(spectrum_int_spatial.value / grad(spectrum_int_spatial.frequency))
    inter_freq = uconvert.(unit(_TF_temporal), spectrum_int_spatial.frequency, dispersion)

    spectrum_int_temporal = DiscreteOmnidirectionalSpectrum(inter_value, inter_freq)
    return convert_frequency(spectrum_int_temporal, TF_new, dispersion)
end

"""
    convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF}, TF_new, dispersion::Dispersion=deepwater_gradient)

Converts the spectra into the new frequency units using the [`DimensionfulAngles.Dispersion`](@extref) relation
and returns a new struct with updated spectrum and frequency.

See also [`DimensionfulAngles.Dispersion`](@extref)

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1.0u"Hz", 5.0u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v,f);

julia> s2 = convert_frequency(s1, 1.0u"Hz^-1");

```
"""
function convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS,TF}, TF_new::T,
    dispersion::Equivalence=deepwater_gradient) where {TS,TF<:_Temporal,T<:_Spatial}
    
    spectrum_int_temporal = convert_frequency(spectrum, _TF_temporal, dispersion)

    grad = dispersion.gradient_inverse
    inter_value = @. upreferred(spectrum_int_temporal.value / grad(spectrum_int_temporal.frequency))
    inter_freq = uconvert.(unit(_TF_spatial), spectrum_int_temporal.frequency, dispersion)

    spectrum_int_spatial = DiscreteOmnidirectionalSpectrum(inter_value, inter_freq)
    return convert_frequency(spectrum_int_spatial, TF_new, dispersion)
end

