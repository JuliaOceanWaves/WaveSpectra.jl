struct WaveTimeSeries{TS<:Quantity, TF<:Quantity, N} <: AbstractArray{TS, N}
    value::AbstractVecOrMat{<:Quantity}
    time::AbstractVector{<:Quantity}
    colnames::Union{AbstractVector{<:AbstractString}, Nothing}
    meta::Any
    function WaveTimeSeries(value::AbstractVecOrMat{<:Quantity}, time::AbstractVector{<:Quantity}, colnames=nothing, meta=nothing)
        # Parameters
        N=ndims(value) # Don't need to check N if AbstractVecOrMat handles 0 dim and 3+ dims
        TF=eltype(time)
        TS=eltype(value)
        # Check Arguments
        (size(value)[1] ≠ length(time)) && error("'time' and first dimension of 'value' array must be same size")
        return new{TS,TF,N}(value, time)
    end
end

#Constructors
function WaveTimeSeries(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF}) where {TS, TF}
    v, f = spectrum.value, spectrum.frequency
    N = size(v, 1)
    tend = 2π/f[1]
    ϕ = -im .* randn(N)
    V = @. exp(ϕ) * √(2*v*(f[2] - f[1]))
    X = N .*irfft(ustrip.(V), 2*N-1)
    t = range(0/unit(TF), tend, 2*N-1)
    return WaveTimeSeries(X.*unit(eltype(V)), t)
end

_firstState() = 1
#Interface
Base.size(series::WaveTimeSeries) = size(series.value)
Base.getindex(series::WaveTimeSeries, i::Int) = series.value[i]
Base.getindex(series::WaveTimeSeries, I::Vararg{Int, 2}) = series.value[I...]
Base.getindex(series::WaveTimeSeries, I) = [series.value[i] for i in I]
Base.length(series::WaveTimeSeries) = length(series.value)
Base.iterate(series::WaveTimeSeries) = isempty(series.value) ? nothing : (series.value[_firstState()], _firstState() + 1)
Base.iterate(series::WaveTimeSeries, state::Int) = (state <= length(series.value)) ? (series.value[state], state+1) : nothing
Base.firstindex(::WaveTimeSeries) = _firstState()
Base.lastindex(series::WaveTimeSeries) = length(series.value)
Base.eltype(::WaveTimeSeries{TS}) where {TS} = TS


function Base.show(io::IO, series::WaveTimeSeries{TS,TF,N}) where {TS,TF,D,N}
    print(io, "\nPARAMS: {$TS,\n\t $TF,\n\t Value Dims=$N}\nTIME {$(unit(TF))}, VALUES $(size(series))\t{$(unit(TS))}")
end

function Base.show(io::IO, ::MIME"text/plain", series::WaveTimeSeries{TS,TF,N}) where {TS,TF,N}
    if typeof(series.time) <: AbstractRange
        af = ustrip.(series.time)
    else
        af = round.(ustrip.(series.time), digits=3)
        af = length(af) > 2 ? join([af[begin], af[end]], " … ") : af
    end

    if N < 2
        av = round.(ustrip.(series.value), digits=3)
        av = length(av) > 2 ? join([av[begin], av[end]], " … ") : av
    else
        lrow_str = "\n\t\t\t⋮   ⋱   ⋮\n\t\t\t"
        erow = [round.(r,digits=3) for r in eachrow(ustrip.(series.value))]
        rows = [length(r) > 2 ? join([r[begin], r[end]], " … ") : join(r, ", ") for r in erow]
        av = length(rows) > 2 ? join([rows[begin], rows[end]], lrow_str) : join(rows, ", \n")
    end
    print(io, "PARAMS: {$TS,\n\t $TF,\n\t Value Dims=$N}\nTIME\t{$(unit(TF))}: \t$(af)\nVALUES $(size(series))\t{$(unit(TS))}: \t$(av)")
end

@recipe function f(series::WaveTimeSeries, args...)
    xlabel --> "Time"
    ylabel --> "Wave Elevation"
    marker := :auto
    (series.time, series.value, args...)
end