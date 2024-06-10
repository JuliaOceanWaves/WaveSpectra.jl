# Structs definition





# Interface
function isunitful(spectrum::OmnidirectionalSpectrum{S, F})  where {S, F}
    return isa(1*upreferred(F()), Frequency)
end

function isunitful(spectrum::DiscreteOmnidirectionalSpectrum{S,F,D,M}) where {S,F,D,M}
    return isa(1*upreferred(F()), Frequency)
end

function Unitful.unit(spectrum::OmnidirectionalSpectrum{S, F})  where {S, F}
    return unit(spectrum(0*upreferred(F())))
end

function Unitful.unit(spectrum::DiscreteOmnidirectionalSpectrum{S,F,D,M})  where {S,F,D,M}
    return unit(eltype(spectrum))
end

Unitful.dimension(spectrum::OmnidirectionalSpectrum{S, F}) where {S, F} = S()

Unitful.dimension(spectrum::DiscreteOmnidirectionalSpectrum{S,F,D,M}) where {S,F,D,M} = S()

function quantity(spectrum::OmnidirectionalSpectrum{S, F}) where {S, F}
    dims = S()*F()
    units = upreferred(dims)
    return dims, units
end

function quantity(spectrum::DiscreteOmnidirectionalSpectrum{S,F,true,M}) where {S,F,M}
    dims = S()*F()
    units = upreferred(dims)
    return dims, units
end

function quantity(spectrum::DiscreteOmnidirectionalSpectrum{S,F,false,M}) where {S,F,M}
    dims = S()
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

Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, i::Int) = spectrum.value[i]

function Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, I::Vararg{Int, ndims(A)})
    return spectrum.value[I...]
end

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

@recipe function f(spectrum::OmnidirectionalSpectrum{S, F, D}, args...) where {S, F, D}
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    if isdiscrete(spectrum)
        marker := :auto
        (spectrum.discrete.frequency, spectrum.discrete.value, args...)
    else
        isempty(args) && (args=(0*upreferred(F()), 10*upreferred(F())))
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

function _update_domain(spectrum::OmnidirectionalSpectrum{S, F, D}, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing) where {S, F, D}
    isnothing(f_begin) && (f_begin = 0*upreferred(F()))
    isnothing(f_end) && (f_end = Inf*upreferred(F()))
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
    spectrum::OmnidirectionalSpectrum{S, F, D},
        f_begin::Union{Number, Nothing}=nothing,
        f_end::Union{Number, Nothing}=nothing;
        alg::AbstractIntegralAlgorithm=QuadGKJL(),
        kwargs...
    ) where {S, F, D}
    f_begin, f_end = _update_domain(spectrum, f_begin, f_end)
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8*upreferred(S()*F())
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
        spectrum::AbstractMatrix{<:Number}, frequency::AbstractVector{<:Number},
        n::Integer, method::AbstractSampledIntegralAlgorithm=TrapezoidalRule(); kwargs...
    )
    return integrate(spectrum.*(frequency.^n)', frequency, method; kwargs...)
end

function spectral_moment(spectrum::OmnidirectionalSpectrum{S, F, true}, n::Int, f_begin::Union{Number, Nothing}=nothing,
        f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL();
        kwargs...
    ) where {S, F}
    f_begin, f_end = _update_domain(spectrum, f_begin, f_end)
    if :abstol ∉ keys(kwargs)
        abstol = 1e-8*upreferred(S()*F()^(n+1))
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

function energy_period(
        spectrum::OmnidirectionalSpectrum{S, F, true},
        f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing,
        alg::AbstractIntegralAlgorithm=QuadGKJL(); kwargs...
    ) where {S, F}
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

function significant_waveheight(spectrum::OmnidirectionalSpectrum{S, F, true}, f_begin::Union{Number, Nothing}=nothing, f_end::Union{Number, Nothing}=nothing, alg::AbstractIntegralAlgorithm=QuadGKJL(); kwargs...)
	    m_0 = spectral_moment(spectrum, 0, f_begin, f_end, alg; kwargs...
    ) where {S, F}
    return 4√m_0
end

function peak_period(spectrum::OmnidirectionalSpectrum{S, F, D})

# parametric
# TODO: JONSWAP, PM, Etc
