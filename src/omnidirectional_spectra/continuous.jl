
# Struct
"""
    OmnidirectionalSpectrum(func::Function, TS::DataType, TF::DataType)
    OmnidirectionalSpectrum(func::Function, TF::DataType=typeof(1.0Hz))
    OmnidirectionalSpectrum(value::AbstractVector{<:Quantity}, frequency::AbstractVector{<:Quantity}; interpolation::Function=linear_interpolation)
    OmnidirectionalSpectrum(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, true, 1}; interpolation::Function=linear_interpolation) where {TS, TF}

Create an omnidirectional spectrum using a combination of a function with the appropriate 
data types, two vectors, or an existing DiscreteOmnidirectionalSpectrum.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s1 = OmnidirectionalSpectrum(func, typeof(1u"Hz"), typeof(1u"Hz"));

julia> s2 = OmnidirectionalSpectrum(func, typeof(1u"Hz"));

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s3 = OmnidirectionalSpectrum(v,f);

julia> sd = DiscreteOmnidirectionalSpectrum(v,f);

julia> s4 = OmnidirectionalSpectrum(sd);

```
"""
struct OmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity} <: Function
    func :: Function

    function OmnidirectionalSpectrum(func::Function, TS::DataType, TF::DataType)
        (dimension(TF) ∉ _frequency_dims) && error("invalid frequency dimensions")
        (typeof(func(ones(TF)[])) ≠ TS) && error("invalid spectrum dimensions")
        return new{TS, TF}(func)
    end
end

# Constructors
function OmnidirectionalSpectrum(func::Function, TF::DataType=typeof(1.0Hz))
    TS = typeof(func(ones(TF)[]))
    return OmnidirectionalSpectrum(func, TS, TF)
end

function OmnidirectionalSpectrum(
        value::AbstractVector{<:Quantity}, frequency::AbstractVector{<:Quantity};
        interpolation::Function=linear_interpolation
    )
    TS = eltype(value)
    TF = eltype(frequency)
    func = x -> interpolation(frequency, value; extrapolation_bc = 0*unit(TS))(x)
    return OmnidirectionalSpectrum(func, TS, TF)
end

function OmnidirectionalSpectrum(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, true, 1}; 
        interpolation::Function=linear_interpolation) where {TS, TF}
    return OmnidirectionalSpectrum(spectrum.value, spectrum.frequency; interpolation)
end

# Call Methods
function (spectrum::OmnidirectionalSpectrum{TS,TF})(frequency::Quantity) where {TS,TF}
    (typeof(frequency) ≠ TF) && (frequency=convert(TF, frequency))
    return spectrum.func(frequency)
end

# Unitful interface
"""
    Unitful.unit(::OmnidirectionalSpectrum{TS, TF})

Extends from Unitful, returns the units of the spectra values.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"))
(::OmnidirectionalSpectrum{Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}, Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}})     (generic function with 1 method)

julia> unit(s)
Hz

```
"""
Unitful.unit(::OmnidirectionalSpectrum{TS, TF})  where {TS, TF} = unit(TS)
"""
    Unitful.unit(::OmnidirectionalSpectrum{TS, TF})

Extends from Unitful, returns the dimension of the spectra values.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"))
(::OmnidirectionalSpectrum{Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}, Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}})     (generic function with 1 method)

julia> dimension(s)
𝐓⁻¹

```
"""
Unitful.dimension(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TS)
"""
    Unitful.unit(::OmnidirectionalSpectrum{TS, TF})

Extends from Unitful, returns the units of the expected frequency.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"))
(::OmnidirectionalSpectrum{Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}, Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}})     (generic function with 1 method)

julia> unit(s)
Hz

```
"""
frequency_unit(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = unit(TF)
"""
    Unitful.unit(::OmnidirectionalSpectrum{TS, TF})

Extends from Unitful, returns the dimension of the expected frequency.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"))
(::OmnidirectionalSpectrum{Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}, Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}})     (generic function with 1 method)

julia> dimension(s)
𝐓⁻¹

```
"""
frequency_dimension(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TF)

"""
    quantity(::OmnidirectionalSpectrum{TS, TF})

Return the dimensions and units of the product between spectra and frequency.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"))
(::OmnidirectionalSpectrum{Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}, Quantity{Float64, 𝐓⁻¹, Unitful.FreeUnits{(Hz,), 𝐓⁻¹, nothing}}})     (generic function with 1 method)

julia> quantity(s)
(𝐓⁻², Hz²)

```
"""
function quantity(::OmnidirectionalSpectrum{TS, TF}) where {TS, TF}
    dimensions = dimension(TS) * dimension(TF)
    units = unit(TS) * unit(TF)
    return dimensions, units
end

# Plots recipe
@recipe function f(spectrum::OmnidirectionalSpectrum{TS, TF}, kwargs...) where {TS, TF}
    xlabel --> "frequency"
    ylabel --> "spectral density"
    return (spectrum.func, 0unit(TF), 10unit(TF), kwargs...)
end

@recipe function f(
        spectrum::OmnidirectionalSpectrum{TS,TF}, xmin::Quantity, xmax::Quantity, kwargs...
    ) where {TS,TF}
    xlabel --> "frequency"
    ylabel --> "spectral density"
    return (spectrum.func, xmin, xmax, kwargs...)
end

# Spectral moments
"""
    spectrum::OmnidirectionalSpectrum{TS, TF}, n::Int, f_begin::Union{Quantity, Nothing}=nothing, 
        f_end::Union{Quantity, Nothing}=nothing; alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...

Calculate the n-th spectral moment of a discrete spectra struct, use given range if specified.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> func(x) = x
func (generic function with 1 method)

julia> s = OmnidirectionalSpectrum(func, typeof(1.0u"Hz"));

julia> spectral_moment(s, -1, 1u"Hz", 5u"Hz")
4.0 s⁻¹

julia> spectral_moment(s, 0, 1u"Hz", 5u"Hz")
12.0 s⁻²

julia> spectral_moment(s, 1, 1u"Hz", 5u"Hz")
41.333333333333336 s⁻³
```
"""
function spectral_moment(spectrum::OmnidirectionalSpectrum{TS, TF}, n::Int,
        f_begin::Union{Quantity, Nothing}=nothing, f_end::Union{Quantity, Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...) where {TS,TF}
    isnothing(f_begin) && (f_begin = 0 * unit(TF))
    isnothing(f_end) && (f_end = Inf * unit(TF))
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8 * unit(TS) * unit(TF)^(n+1)
        kwargs = merge(values(kwargs), (abstol=abstol,))
    end
    sol = solve(IntegralProblem((f, _) -> spectrum(f)*f^n, (f_begin, f_end), nothing),
        alg; kwargs...)
    if sol.retcode ≠ ReturnCode.Success
        error("solution unsuccessful with code: $(sol.retcode)")
    end
    return upreferred(sol.u)
end

# Convert frequency
function convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new,
    dispersion::Dispersion=Dispersion()) where {TS,TF}
    grad = _grad[(dimension(TF), dimension(TF_new))]
    function func(f)
        f_org = uconvert(unit(TF), f, dispersion)
        return upreferred(spectrum.func(f_org) / grad(f_org))
    end
    return OmnidirectionalSpectrum(func, TF_new)
end
