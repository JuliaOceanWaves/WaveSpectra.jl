# TODO: Display of Spectrum and OmnidirectionalSpectrum in Pluto... combination of DataFrames and AxisArrays
# TODO: Option to display as images in html?

# 2D Spectrum
function Base.summary(io::IO, z::Spectrum)
    println(io, "Spectrum of quantity in [", unit(z, :integral), "] with axes:")
    println(io, "    ", String(axistype(z.axis1)), ", ", Vector{Quantity{eltype(z.axis1).parameters[1]}}(z.axis1))
    println(io, "    ", String(axistype(z.axis2)),  ", ", Vector{Quantity{eltype(z.axis2).parameters[1]}}(z.axis2))
    print(io, "And data, a ", join(string.(size(z)), "x"), " Matrix{Quantity{", eltype(z).parameters[1],"}}")
    return nothing
end

# arrays of 2D Spectra
function Base.summary(io::IO, z::Array{<:Spectrum})
    print(io, join(string.(size(z)),"x"), " Array{Spectrum, ", string(ndims(z)), "}")
end

function Base.summary(io::IO, z::Matrix{<:Spectrum})
    print(io, join(string.(size(z)),"x"), " Matrix{Spectrum}")
end

function Base.summary(io::IO, z::Vector{<:Spectrum})
    print(io, string.(length(z)), "-element Vector{Spectrum}")
end

# AxisArray: array of 2D Spectra
function Base.summary(io::IO, A::AxisArray{<:Spectrum,N}) where N
    println(io, "$N-dimensional AxisArray{Spectrum,$N,…} with axes:")
    for (name, val) in zip(axisnames(A), axisvalues(A))
        print(io, "    :$name, ")
        show(IOContext(io, :limit=>true), val)
        println(io)
    end
    print(io, "And data, a ", summary(A.data))
end

# AxisArray: 2D Spectrum
function Base.summary(io::IO, A::AxisArray{T, N, <:Spectrum}) where {T,N}
    println(io, "$N-dimensional AxisArray{Quantity{…}, $N, Spectrum{…}, …} with axes:")
    for (name, val) in zip(axisnames(A), axisvalues(A))
        println(io, "    :$(name),  $(Vector{Quantity{eltype(val).parameters[1]}}(val))")
    end
    print(io, "And data, a Spectrum")
end


# 1D
function Base.summary(io::IO, z::OmnidirectionalSpectrum)
    _summary = """
        OmnidirectionalSpectrum of quantity in [$(unit(z, :integral))] with axis:
            $(String(axistype(z.axis))), $(Vector{Quantity{eltype(z.axis).parameters[1]}}(z.axis))
        And data, a $(length(z.axis))-element Vector{Quantity{$(eltype(z).parameters[1])}}"""
    print(io, _summary)
end

# arrays of 1D spectra
function Base.summary(io::IO, z::Array{<:OmnidirectionalSpectrum})
    print(io, join(string.(size(z)),"x"), " Array{OmnidirectionalSpectrum, ", string(ndims(z)), "}")
end

function Base.summary(io::IO, z::Matrix{<:OmnidirectionalSpectrum})
    print(io, join(string.(size(z)),"x"), " Matrix{OmnidirectionalSpectrum}")
end

function Base.summary(io::IO, z::Vector{<:OmnidirectionalSpectrum})
    print(io, string.(length(z)), "-element Vector{OmnidirectionalSpectrum}")
end

# AxisArray: array of 1D Spectra
function Base.summary(io::IO, A::AxisArray{<:OmnidirectionalSpectrum,N}) where N
    println(io, "$N-dimensional AxisArray{OmnidirectionalSpectrum,$N,…} with axis:")
    for (name, val) in zip(axisnames(A), axisvalues(A))
        print(io, "    :$name, ")
        show(IOContext(io, :limit=>true), val)
        println(io)
    end
    print(io, "And data, a ", summary(A.data))
end

# AxisArray: 1D Spectrum
function Base.summary(io::IO, A::AxisArray{T, N, <:OmnidirectionalSpectrum}) where {T,N}
    println(io, "$N-dimensional AxisArray{Quantity{…}, $N, OmnidirectionalSpectrum{…}, …} with axis:")
    name, val = axisnames(A)[1], axisvalues(A)[1]
    println(io, "    :$(name),  $(Vector{Quantity{eltype(val).parameters[1]}}(val))")
    print(io, "And data, an OmnidirectionalSpectrum")
end
