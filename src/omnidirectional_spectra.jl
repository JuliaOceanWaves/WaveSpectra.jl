# Struct definition
struct OmnidirectionalSpectrum{S, F, D}
    func::Function
    discrete::Union{Nothing, @NamedTuple begin
        frequency::AbstractVector
        value::AbstractVector
    end}

    function OmnidirectionalSpectrum{S, F, D}(func, discrete) where {S<:Dimensions, F<:Dimensions, D}
        # check parameters
        !isa(D, Bool) && error("parameter `D` must be a boolean")
        if !(isa(1*upreferred(F()), Frequency) || (F == typeof(NoDims)))
            error("parameter `F` must be a `Unitful.Frequency` or `Unitful.NoDims`")
        end
        if (F == typeof(NoDims)) && !(S == typeof(NoDims))
            error("if `F` is dimensionless `S` must be too")
        end
        # check arguments
        if !(isnothing(discrete))
            if (length(discrete.frequency) !== length(discrete.value))
                error("`frequency` and `value` arrays must have same size")
            end
            if func.(discrete.frequency) ≉ discrete.value
                error("provided function and values do not agree")
            end
            ftype = eltype(discrete.frequency)
            stype = eltype(discrete.value)
            if !(ftype<:Frequency || dimension(ftype)==NoDims)
                error("`discrete.frequency` must have a frequency dimension or no dimension")
            end
            if dimension(ftype)==NoDims && !(dimension(stype)==NoDims)
                error("if `discrete.frequency` is dimensionless, `discrete.value` must be too")
            end
        end
        # create parametric object
        return new{S, F, D}(func, discrete)
    end
end

# Constructors
function OmnidirectionalSpectrum(func::Function; dims_frequency::Dimensions=𝐓^-1, dims_value::Dimensions=𝐋^2*𝐓, density=true)
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, nothing)
end

function OmnidirectionalSpectrum(func::Function, freqs::AbstractVector{<:Number}; density=true)
    dims_frequency = dimension(eltype(freqs))
    dims_value = dimension(func(0*upreferred(eltype(freqs))))
    discrete = (
        frequency = freqs,
        value = func.(frequency)
    )
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, discrete)
end

function OmnidirectionalSpectrum(spectrum::OmnidirectionalSpectrum, freqs::AbstractVector{<:Number})
    Ts, Tf, density = typeof(spectrum).parameters
    discrete = (
            frequency = freqs,
            value = spectrum.func.(frequency)
        )
    return OmnidirectionalSpectrum{Ts, Tf, density}(spectrum.func, discrete)
end

function OmnidirectionalSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; density=true, interpolation::Function=linear_interpolation)
    dims_frequency = dimension(eltype(freqs))
    dims_value = dimension(eltype(vals))
    discrete = (
        frequency = freqs,
        value = vals
    )
    func = x->interpolation(freqs, vals; extrapolation_bc = 0*unit(eltype(vals)))(x)
    @assert func.(freqs) ≈ vals
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, discrete)
end

# Interface
"""
    isdiscrete(spectrum::OmnidirectionalSpectrum)

Return whether the OmnidirectionalSpectrum uses discrete vectors.

# Examples
```jldoctest
a = 1
b = 2
a + b

# output

3
```
"""
isdiscrete(spectrum::OmnidirectionalSpectrum) = !isnothing(spectrum.discrete)
isdensity(spectrum::OmnidirectionalSpectrum) = typeof(spectrum) <: OmnidirectionalSpectrum{S, F, true} where {S,F}
isunitful(spectrum::OmnidirectionalSpectrum{Ts, Tf, D})  where {Ts, Tf, D} = isa(1*upreferred(Tf()), Frequency)
Unitful.unit(spectrum::OmnidirectionalSpectrum{Ts, Tf, D})  where {Ts, Tf, D} = isdiscrete(spectrum) ? unit(eltype(spectrum.discrete.value)) : unit(spectrum(0*upreferred(Tf())))
Unitful.dimension(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}) where {Ts, Tf, D} = Ts()
function quantity(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}) where {Ts, Tf, D}
    dims = Ts()*Tf()
    units = upreferred(dims)
    return dims, units
end

(spectrum::OmnidirectionalSpectrum)(freq::Number) =  spectrum.func(freq)

mutable struct DiscreteError <: Exception
    message::String
    DiscreteError() = new("continuous `OmnidirectionalSpectrum` instance was treated as discrete")
end
Base.showerror(io::IO, e::DiscreteError) = print(io, "DiscreteError: $(e.message)")

_firststate() = 1
Base.getindex(spectrum::OmnidirectionalSpectrum, i::Int) = isdiscrete(spectrum) ? spectrum.discrete.value[i] : throw(DiscreteError())
Base.getindex(spectrum::OmnidirectionalSpectrum, I) = isdiscrete(spectrum) ? [spectrum.discrete.value[i] for i in I] : throw(DiscreteError())
Base.firstindex(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? _firststate() : throw(DiscreteError())
Base.lastindex(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? length(spectrum.discrete.value) : throw(DiscreteError())
Base.length(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? length(spectrum.discrete.value) : throw(DiscreteError())
Base.size(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? size(spectrum.discrete.value) : throw(DiscreteError())
Base.eltype(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? eltype(spectrum.discrete.value) : throw(DiscreteError())
Base.iterate(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? (spectrum.discrete.value[firstindex(spectrum)], _firststate()+1) : throw(DiscreteError())
function Base.iterate(spectrum::OmnidirectionalSpectrum, state::Int)
    isdiscrete(spectrum) || throw(DiscreteError())
    return (state<=length(spectrum)) ? (spectrum.discrete.value[state], state+1) : nothing
end

# Plots recipes
function _labels(spectrum::OmnidirectionalSpectrum)
    x_label = "frequency"
    y_label = isdensity(spectrum) ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}, args...) where {Ts, Tf, D}
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    if isdiscrete(spectrum)
        marker := :auto
        (spectrum.discrete.frequency, spectrum.discrete.value, args...)
    else
        isempty(args) && (args=(0*upreferred(Tf()), 10*upreferred(Tf())))
        (spectrum.func, args...)
    end
end

@recipe function f(spectrum::OmnidirectionalSpectrum, freq::AbstractVector{<:Number}, args...)
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (freq, spectrum.(freq), args...)
end

# integrals/statistics
function _integrate(func::Function, f_begin::Union{Number, Nothing}=nothing,
        f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL();
        kwargs...
    )
    sol = solve(IntegralProblem((x,p)->func(x), (f_begin, f_end), nothing), alg; kwargs...)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return upreferred(sol.u)
end

function _update_domain(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing) where {Ts, Tf, D}
    isnothing(f_begin) && (f_begin = 0*upreferred(Tf()))
    isnothing(f_end) && (f_end = Inf*upreferred(Tf()))
    return (f_begin, f_end)
end

function integrate(y::AbstractVector{<:Number}, x::AbstractVector{<:Number},
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    sol = solve(SampledIntegralProblem(y, x), method; kwargs...)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return upreferred(sol.u)
end

function integrate(y::AbstractMatrix{<:Number}, x::AbstractVector{<:Number},
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    return [integrate(i_y, x, method; kwargs...) for i_y in eachrow(y)]
end

function integrate(
    spectrum::OmnidirectionalSpectrum{Ts, Tf, D},
        f_begin::Union{Number, Nothing}=nothing,
        f_end::Union{Number, Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(),
        kwargs...
    ) where {Ts, Tf, D}
    f_begin, f_end = _update_domain(spectrum, f_begin, f_end)
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8*upreferred(Ts()*Tf())
        kwargs = merge(values(kwargs), (abstol=abstol,))
    end
    return _integrate(spectrum.func, f_begin, f_end, alg; kwargs...)
end

function spectral_moment(
        spectrum::AbstractVector{<:Number}, frequency::AbstractVector{<:Number}, n::Integer,
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    return integrate(spectrum.*frequency.^n, frequency, method; kwargs...)
end

function spectral_moment(
    spectrum::AbstractMatrix{<:Number}, frequency::AbstractVector{<:Number}, n::Integer,
    method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
)
return integrate(spectrum.*(frequency.^n)', frequency, method; kwargs...)
end

function spectral_moment(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}, n::Int, f_begin::Union{Number, Nothing}=nothing,
        f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL();
        kwargs...
    ) where {Ts, Tf, D}
    f_begin, f_end = _update_domain(spectrum, f_begin, f_end)
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8*upreferred(Ts()*Tf()^(n+1))
        kwargs = merge(values(kwargs), (abstol=abstol,))
    end
    return _integrate(f -> spectrum(f)*f^n, f_begin, f_end, alg; kwargs...)
end

function energy_period(
        spectrum::AbstractVecOrMat{<:Number}, frequency::AbstractVector{<:Number},
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    m_n1 = spectral_moment(spectrum, frequency, -1, method; kwargs...)
    m_0 = spectral_moment(spectrum, frequency, 0, method; kwargs...)
    return m_n1./m_0
end

function energy_period(spectrum::OmnidirectionalSpectrum, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL(); kwargs...)
    m_n1 = spectral_moment(spectrum, -1, f_begin, f_end, alg; kwargs...)
    m_0 = spectral_moment(spectrum, 0, f_begin, f_end, alg; kwargs...)
    return m_n1/m_0
end

function significant_waveheight(
        spectrum::AbstractVecOrMat{<:Number}, frequency::AbstractVector{<:Number},
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    m_0 = spectral_moment(spectrum, frequency, 0, method; kwargs...)
    return 4*.√m_0
end

function significant_waveheight(spectrum::OmnidirectionalSpectrum, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL(); kwargs...)
	m_0 = spectral_moment(spectrum, 0, f_begin, f_end, alg; kwargs...)
    return 4√m_0
end

# parametric
# TODO: JONSWAP, PM, Etc
