# Struct
"""
    DiscreteOmnidirectionalSpectrum(value::AbstractVecOrMat{<:Quantity}, frequency::AbstractVector{<:Quantity}; density::Bool=true)
    DiscreteOmnidirectionalSpectrum(func::Function, frequency::AbstractVector{<:Quantity}; density::Bool=true)

Create discrete omnidirectional spectrum Struct using two Unitful vectors, matrix and vector, or function and vector.

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

julia> func(x) = x
func (generic function with 1 method)

julia> s3 = DiscreteOmnidirectionalSpectrum(func, f);

```
"""
struct DiscreteOmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity, D, N} <: AbstractArray{TS, N}
    value::AbstractVecOrMat{<:Quantity}
    frequency::AbstractVector{<:Quantity}
    function DiscreteOmnidirectionalSpectrum(value::AbstractVecOrMat{<:Quantity}, frequency::AbstractVector{<:Quantity}; density::Bool=true)
        # Parameters
        N=ndims(value)
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

# Call Methods
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

using Unitful:ustrip
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
"""
    Unitful.unit(::DiscreteOmnidirectionalSpectrum)

Extends from Unitful, returns the units of the spectra vector/matrix in the struct

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=range(1u"m/Hz", 3u"m/Hz", 3)
(1.0:1.0:3.0) m Hz⁻¹

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> unit(s)
m Hz⁻¹

```
"""
Unitful.unit(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = unit(TS)
"""
    Unitful.dimension(::DiscreteOmnidirectionalSpectrum)

Extends from Unitful, returns the dimension of the spectra vector/matrix in the struct

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=range(1u"m/Hz", 3u"m/Hz", 3)
(1.0:1.0:3.0) m Hz⁻¹

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> dimension(s)
𝐋 𝐓

```
"""
Unitful.dimension(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = dimension(TS)
"""
    frequency_unit(::DiscreteOmnidirectionalSpectrum)

Return the units of the frequency vector in the struct

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=range(1u"m/Hz", 3u"m/Hz", 3)
(1.0:1.0:3.0) m Hz⁻¹

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> frequency_unit(s)
Hz

```
"""
frequency_unit(::DiscreteOmnidirectionalSpectrum{TS, TF}) where {TS, TF} = unit(TF)
"""
    frequency_dimension(::DiscreteOmnidirectionalSpectrum)

Return the dimension of the frequency vector in the struct

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=range(1u"m/Hz", 3u"m/Hz", 3)
(1.0:1.0:3.0) m Hz⁻¹

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> frequency_dimension(s)
𝐓⁻¹

```
"""
frequency_dimension(::DiscreteOmnidirectionalSpectrum{TS, TF}) where {TS, TF} = dimension(TF)

"""
    quantity(::DiscreteOmnidirectionalSpectrum)

Return the dimensions and units of the product between spectra and frequency.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=range(1u"m/Hz", 3u"m/Hz", 3)
(1.0:1.0:3.0) m Hz⁻¹

julia> f=range(1u"Hz", 3u"Hz", 3)
(1.0:1.0:3.0) Hz

julia> s = DiscreteOmnidirectionalSpectrum(v,f);

julia> quantity(s)
(𝐋, m)

```
"""
function quantity(::DiscreteOmnidirectionalSpectrum{TS, TF, D}) where {TS, TF, D}
    dimensions = dimension(TS) * dimension(TF)
    units = unit(TS) * unit(TF)
    return dimensions, units
end

# Plots recipes
function _labels(::DiscreteOmnidirectionalSpectrum{TS, TF, D}) where {TS, TF, D}
    x_label = "frequency"
    y_label = D ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(spectrum::DiscreteOmnidirectionalSpectrum, args...)
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (spectrum.frequency, spectrum.value, args...)
end

"""
    spectral_moment(spectrum::DiscreteOmnidirectionalSpectrum, n::Real=0; alg::AbstractIntegralAlgorithm=TrapezoidalRule())

Calculate the n-th spectral moment of a discrete spectra struct.

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

"""
    convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF}, TF_new, dispersion::Dispersion=Dispersion()) where {TS, TF}

Converts the spectra into the new frequency units using the [`DimensionfulAngles.Dispersion`](@extref) relation 
and returns a new struct with updated spectrum and frequency.

See also [`DimensionfulAngles.Dispersion`](@extref)

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s1 = DiscreteOmnidirectionalSpectrum(v,f);

julia> s2 = convert_frequency(s1, typeof(1u"Hz^-1"));

```
"""
function convert_frequency(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF}, TF_new,
    dispersion::Dispersion=Dispersion()) where {TS, TF}
    grad = _grad[(dimension(TF), dimension(TF_new))]
    new_value = @. upreferred(spectrum.value / grad(spectrum.frequency))
    new_freq = uconvert.(unit(TF_new), spectrum.frequency, dispersion)
    return DiscreteOmnidirectionalSpectrum(new_value, new_freq)
end
