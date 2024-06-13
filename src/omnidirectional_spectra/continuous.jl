
# Struct
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

# Call Methods
function (spectrum::OmnidirectionalSpectrum{TS,TF})(frequency::Quantity) where {TS,TF}
    (typeof(frequency) ≠ TF) && (frequency=convert(TF, frequency))
    return spectrum.func(frequency)
end

# Unitful interface
Unitful.unit(::OmnidirectionalSpectrum{TS, TF})  where {TS, TF} = unit(TS)
Unitful.dimension(::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TS)
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
function phase_velocity(frequency::Quantity, dispersion::Dispersion)
    wavelength = uconvert(m, frequency, dispersion)
    period = uconvert(s, frequency, dispersion)
    return wavelength/period
end

function convert_frequency(spectrum::OmnidirectionalSpectrum{TS,TF}, TF_new,
    dispersion::Dispersion=Dispersion()) where {TS,TF}
    # grad = _get_grad(dimension(TF), dimension(TF_new), dispersion)
    f -> phase_velocity(f, dispersion)



    function func(f)
        f_org = uconvert(unit(TF), f, dispersion)
        return upreferred(spectrum.func(f_org) / grad(f_org))
    end
    return OmnidirectionalSpectrum(func, TF_new)
end
