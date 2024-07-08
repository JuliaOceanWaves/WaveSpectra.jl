struct WaveTimeSeries{TS<:Quantity, TF<:Quantity, N} <: AbstractArray{TS, N}
    time::AbstractVector{<:Quantity}
    values::AbstractVecOrMat{<:Quantity}
    colnames::Union{AbstractVector{<:AbstractString}, Nothing}
    meta::Any
    function WaveTimeSeries(values::AbstractVecOrMat{<:Quantity}, time::AbstractVector{<:Quantity}, colnames=nothing, meta=nothing)
        # Parameters
        N=ndims(values) # Don't need to check N if AbstractVecOrMat handles 0 dim and 3+ dims
        TF=eltype(time)
        TS=eltype(values)
        # Check Arguments
        (size(values)[1] ≠ length(time)) && error("'time' and first dimension of 'values' array must be same size")
        return new{TS,TF,N}(values, time)
    end
end

#Constructors
function WaveTimeSeries(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, D, 1}; r_amps=nothing, seed=nothing) where {TS, TF, D}
    values, freqs = spectrum.values, spectrum.frequency
    # @assert typeof(freqs) <: AbstractRange "Frequencies must be a range, unevenly spaces frequencies are not yet implemented!"
    n_samples = size(values, 1)
    tend = isapprox(0.0, freqs[begin]) ? 2π/freqs[2] : 2π/freqs[begin]
    t = range(0/unit(TF), tend, 2*n_samples-1)
    if isnothing(seed)
        ϕ = -im .* randn(n_samples)
    else
        ϕ = -im .* randn(seed, n_samples)
    end
    Δfreq = typeof(freqs) <: AbstractRange ? step(freqs) : (isapprox(0.0, freqs[begin]) ? freqs[2] : freqs[begin])
    V = @. exp(ϕ) * sqrt(2*values*Δfreq)
    X = n_samples .*irfft(ustrip.(V), 2*n_samples-1)
    return WaveTimeSeries(X.*unit(eltype(V)), t)
end

#Interface
Base.size(series::WaveTimeSeries) = size(series.values)
Base.getindex(series::WaveTimeSeries, i::Int) = series.values[i]
Base.getindex(series::WaveTimeSeries, I::Vararg{Int, 2}) = series.values[I...]
Base.getindex(series::WaveTimeSeries, I) = [series.values[i] for i in I]
Base.length(series::WaveTimeSeries) = length(series.values)
Base.iterate(series::WaveTimeSeries) = isempty(series.values) ? nothing : (series.values[_firstState()], _firstState() + 1)
Base.iterate(series::WaveTimeSeries, state::Int) = (state <= length(series.values)) ? (series.values[state], state+1) : nothing
Base.firstindex(::WaveTimeSeries) = _firstState()
Base.lastindex(series::WaveTimeSeries) = length(series.values)
Base.eltype(::WaveTimeSeries{TS}) where {TS} = TS


function Base.show(io::IO, series::WaveTimeSeries{TS,TF,N}) where {TS,TF,N}
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
        av = round.(ustrip.(series.values), digits=3)
        av = length(av) > 2 ? join([av[begin], av[end]], " … ") : av
    else
        lrow_str = "\n\t\t\t⋮   ⋱   ⋮\n\t\t\t"
        erow = [round.(r,digits=3) for r in eachrow(ustrip.(series.values))]
        rows = [length(r) > 2 ? join([r[begin], r[end]], " … ") : join(r, ", ") for r in erow]
        av = length(rows) > 2 ? join([rows[begin], rows[end]], lrow_str) : join(rows, ", \n")
    end
    print(io, "PARAMS: {$TS,\n\t $TF,\n\t Value Dims=$N}\nTIME\t{$(unit(TF))}: \t$(af)\nVALUES $(size(series))\t{$(unit(TS))}: \t$(av)")
end

@recipe function f(series::WaveTimeSeries, args...)
    xlabel --> "Time"
    ylabel --> "Wave Elevation"
    marker := :auto
    (series.time, series.values, args...)
end