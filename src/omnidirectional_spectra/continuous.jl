# Struct
"""
    OmnidirectionalSpectrum(func::Function, TS::DataType, TF::DataType)

Create an omnidirectional spectrum using a a function with the appropriate data types.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> s1 = OmnidirectionalSpectrum(x -> x, typeof(1.0u"Hz"), typeof(1.0u"Hz"));
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
"""
    OmnidirectionalSpectrum(func::Function, TF::DataType=typeof(1.0Hz))

Create an omnidirectional spectrum using a function with the appropriate input data type.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> s = OmnidirectionalSpectrum(x -> x, typeof(1.0u"Hz"));

```
"""
function OmnidirectionalSpectrum(func::Function, TF::DataType=typeof(1.0Hz))
    TS = typeof(func(ones(TF)[]))
    return OmnidirectionalSpectrum(func, TS, TF)
end

"""
    OmnidirectionalSpectrum(value::AbstractVector{<:Quantity}, frequency::AbstractVector{<:Quantity}; interpolation::Function=linear_interpolation)

Create an omnidirectional spectrum using a two vectors.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> s = OmnidirectionalSpectrum(v,f);

```
"""
function OmnidirectionalSpectrum(
        value::AbstractVector{<:Quantity}, frequency::AbstractVector{<:Quantity};
        interpolation::Function=linear_interpolation
    )
    TS = eltype(value)
    TF = eltype(frequency)
    func = x -> interpolation(frequency, value; extrapolation_bc = 0*unit(TS))(x)
    return OmnidirectionalSpectrum(func, TS, TF)
end

"""
    OmnidirectionalSpectrum(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, true, 1}; interpolation::Function=linear_interpolation)

Create an omnidirectional spectrum using an existing [`DiscreteOmnidirectionalSpectrum`](@ref).

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> v=f=range(1u"Hz", 5u"Hz", 5)
(1.0:1.0:5.0) Hz

julia> sd = DiscreteOmnidirectionalSpectrum(v,f);

julia> s = OmnidirectionalSpectrum(sd);

```
"""
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
Unitful.unit(::OmnidirectionalSpectrum{TS, TF})  where {TS, TF} = unit(TS)
Unitful.dimension(::OmnidirectionalSpectrum{TS,TF}) where {TS, TF} = dimension(TS)
frequency_unit(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = unit(TF)
frequency_dimension(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TF)
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
    spectral_moment(spectrum::OmnidirectionalSpectrum{TS, TF}, n::Int,
    f_begin::Union{Quantity, Nothing}=nothing, f_end::Union{Quantity, Nothing}=nothing;
    alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...)
Calculate the n-th spectral moment of a spectra, use range if given.

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> s = OmnidirectionalSpectrum(x -> x, typeof(1.0u"Hz"));

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
        alg::AbstractIntegralAlgorithm=QuadGKJL(), flip::Bool=false, kwargs...) where {TS,TF}
    isnothing(f_begin) && (f_begin = 0 * unit(TF))
    isnothing(f_end) && (f_end = Inf * unit(TF))
    flip && ((f_begin, f_end) = (f_end, f_begin))
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8 * unit(TS) * unit(TF)^(n+1)
        kwargs = merge(values(kwargs), (abstol=abstol,))
    end
    sol = solve(IntegralProblem((f, _) -> spectrum.(f).*f.^n, (f_begin, f_end), nothing),
        alg; kwargs...)
    if sol.retcode ≠ ReturnCode.Success
        error("solution unsuccessful with code: $(sol.retcode)")
    end
    return upreferred(sol.u)
end

for T1 in [_Temporal, _Spatial]
    @eval begin
        function convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new::T,
                dispersion::Equivalence=deepwater_gradient) where {TS,TF<:$T1, T<:$T1}
            # println("Both Temporal/Spatial")
            (1/upreferred(unit(T)) == 1*upreferred(unit(TF))) && @warn "Conversion contains a reciprocal, please use flip keyword in other functions to avoid errors!"
            grad = _get_grad(dimension(TF), dimension(T))
            function func(f)
                f_org = uconvert(unit(TF), f, dispersion)
                return upreferred(spectrum.func(f_org) / grad(f_org))
            end
            return OmnidirectionalSpectrum(func, T)
        end
    end
end

function convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new::T,
            dispersion::Equivalence=deepwater_gradient) where {TS,TF<:_Spatial,T<:_Temporal}
    (1/upreferred(unit(T)) == 1*upreferred(unit(TF))) && @warn "Conversion contains a reciprocal, please use flip keyword in other functions to avoid errors!"
    spectrum_int_spatial = convert_frequency(spectrum, _TF_spatial, dispersion)

    grad = dispersion.gradient
    function func(f)
        f_org = uconvert(unit(_TF_spatial), f, dispersion)
        return upreferred(spectrum_int_spatial.func(f_org) / grad(f_org))
    end

    spectrum_int_temporal = OmnidirectionalSpectrum(func, typeof(_TF_temporal))
    return convert_frequency(spectrum_int_temporal, TF_new, dispersion)
end

"""
    convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new, dispersion::Dispersion=deepwater_gradient)

Converts the spectra into the new frequency units using the [`DimensionfulAngles.Dispersion`](@extref) relation
and returns a new struct with an updated function.

See also [`DimensionfulAngles.Dispersion`](@extref)

# Example
```jldoctest
julia> using WaveSpectra, Unitful

julia> using DimensionfulAngles:radᵃ as rad

julia> s1 = OmnidirectionalSpectrum(x -> x, typeof(1.0u"m/rad"));

julia> s2 = convert_frequency(s1, 1.0u"m");

```
"""
function convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new::T,
            dispersion::Equivalence=deepwater_gradient) where {TS,TF<:_Temporal,T<:_Spatial}
    (1/upreferred(unit(T)) == 1*upreferred(unit(TF))) && @warn "Conversion contains a reciprocal, please use flip keyword in other functions to avoid errors!"
    spectrum_int_temporal = convert_frequency(spectrum, _TF_temporal, dispersion)

    grad = dispersion.gradient_inverse
    function func(f)
        f_org = uconvert(unit(_TF_temporal), f, dispersion)
        return upreferred(spectrum_int_temporal.func(f_org) / grad(f_org))
    end

    spectrum_int_spatial = OmnidirectionalSpectrum(func, typeof(_TF_spatial))
    return convert_frequency(spectrum_int_spatial, TF_new, dispersion)
end
