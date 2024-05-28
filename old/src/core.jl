# TODO: PR for AxisArray to have AbstractAxisArray
# Make Spectrum and OmniSpectrum <: AbstractAxisArray
# Fix problem with AxisArray assignment with Interval

# 2D "full" spectrum

struct Spectrum{TDAT, TAX1, TAX2}  <: AbstractMatrix{TDAT}
    data :: Matrix{TDAT}
    axis1 :: Vector{TAX1}
    axis2 :: Vector{TAX2}

    function Spectrum(data::AbstractMatrix{<:Quantity}, axis1::AbstractVector{<:Quantity}, axis2::AbstractVector{<:Quantity})
        # checks
        @assert size(data) == (length(axis1), length(axis2)) "Data and axes sizes do not match!"
        _check_typeconsistency(data)
        _check_typeconsistency(axis1)
        _check_typeconsistency(axis2)
        (issorted(axis1) && issorted(axis2)) || error("Axes must be increasing.")

        typeaxis1 = axistype(axistype(axis1))[1][1]
        typeaxis2 = axistype(axistype(axis1))[1][1]
        @assert Set((typeaxis1, typeaxis2)) ≠ Set((:direction, :direction))

        # promotion element type
        Tdata = eltype(data).parameters[1]
        Tax1 = eltype(axis1).parameters[1]
        Tax2 = eltype(axis2).parameters[1]
        T = promote_type(Tdata, Tax1, Tax2)

        # convert
        data = Matrix(data .|> T)
        axis1 = Vector(axis1 .|> T)
        axis2 = Vector(axis2 .|> T)

        return new{eltype(data), eltype(axis1), eltype(axis2)}(data, axis1, axis2)
    end
end

function AxisArray(A::Spectrum)
    name1 = axistype(A.axis1)
    name2 = axistype(A.axis2)
    if name1==name2
        name1 = symbol(string(:name1) * "_1")
        name2 = symbol(string(:name2) * "_2")
    end
    axis1 = Axis{name1}(A.axis1)
    axis2 = Axis{name2}(A.axis2)
    return AxisArray(A, axis1, axis2)
end

function Spectrum(A::AxisArray)
    axis1, axis2 = axisvalues(A)
    return Spectrum(A.data, axis1, axis2)
end

Base.size(A::Spectrum{US, U1, U2} where {US, U1, U2}) = size(A.data)
Base.eltype(A::Spectrum{US, U1, U2} where {US, U1, U2}) = eltype(A.data)
Base.copy(A::Spectrum{US, U1, U2} where {US, U1, U2}) = Spectrum(copy(A.data), copy(A.axis1), copy(A.axis2))

Base.getindex(A::Spectrum{US, U1, U2} where {US, U1, U2}, i::Int) = getindex(A.data, i)
Base.getindex(A::Spectrum{US, U1, U2} where {US, U1, U2}, I::Vararg{Int, 2}) = getindex(A.data, I...)
Base.setindex!(A::Spectrum{US, U1, U2} where {US, U1, U2}, v, i::Int) = setindex!(A.data, v, i)
Base.setindex!(A::Spectrum{US, U1, U2} where {US, U1, U2}, v, I::Vararg{Int, 2}) = (A.data[I...] = v)

function Base.getindex(A::Spectrum{US, U1, U2} where {US, U1, U2}; kwargs...)
    kwargs_new = Dict()
    for i in 1:length(kwargs)
        k = keys(kwargs)[i]
        v = kwargs[i]
        kwargs_new[k] = (typeof(v) <: Union{AbstractVector, AbstractInterval}) ? v :
                         typeof(v) <: Integer ? (v:v) :
                         typeof(v) <: Number ? (v .. v) :
                         error("Bad index.")
    end
    return Spectrum(getindex(AxisArray(A); kwargs_new...))
end

function Base.setindex!(A::Spectrum{US, U1, U2} where {US, U1, U2}, v; kwargs...)
    B = AxisArray(A)
    size(v) == (1,1) && (v = v[1, 1]) # scalars
    setindex!(B, v; kwargs...)
    A = Spectrum(B)
    return nothing
end


# 1D "omni" spectrum

struct OmniSpectrum{TDAT, TAX}  <: AbstractVector{TDAT}
    data :: Vector{TDAT}
    axis :: Vector{TAX}

    function OmniSpectrum(data::AbstractVector{<:Quantity}, axis::AbstractVector{<:Quantity})
        # checks
        @assert length(data) == length(axis) "Data and axis lengths do not match!"
        _check_typeconsistency(data)
        _check_typeconsistency(axis)
        issorted(axis) || error("Axis must be increasing.")
        _ = axistype(axis)

        # promotion element type
        Tdata = eltype(data).parameters[1]
        Tax = eltype(axis).parameters[1]
        T = promote_type(Tdata, Tax)

        # convert
        data = Vector(data .|> T)
        axis = Vector(axis .|> T)

        return new{eltype(data), eltype(axis)}(data, axis)
    end
end

function AxisArray(A::OmniSpectrum)
    name = axistype(A.axis)
    axis = Axis{name}(A.axis)
    return AxisArray(A, axis)
end

function OmniSpectrum(A::AxisArray)
    axis = axisvalues(A)
    length(axis)>1 && error("More than one index given for 'OmniSpectrum'.")
    return OmniSpectrum(A.data, axis[1])
end

Base.size(A::OmniSpectrum{US, UA} where {US, UA}) = size(A.data)
Base.eltype(A::OmniSpectrum{US, UA} where {US, UA}) = eltype(A.data)
Base.copy(A::OmniSpectrum{US, UA} where {US, UA}) = OmniSpectrum(copy(A.data), copy(A.axis))

Base.getindex(A::OmniSpectrum{US, UA} where {US, UA}, i::Int) = getindex(A.data, i)
# Base.getindex(A::OmniSpectrum{US, UA} where {US, UA}, I::Vararg{Int, 1}) = getindex(A.data, I...)
Base.setindex!(A::OmniSpectrum{US, UA} where {US, UA}, v, i::Int) = setindex!(A.data, v, i)
# Base.setindex!(A::OmniSpectrum{US, UA} where {US, UA}, v, I::Vararg{Int, 1}) = (A.data[I...] = v)

function Base.getindex(A::OmniSpectrum{US, UA} where {US, UA}; kwargs...)
    length(kwargs)>1 && error("More than one index given for 'OmniSpectrum'.")
    kwargs_new = Dict()
    k = keys(kwargs)[1]
    v = kwargs[1]
    kwargs_new[k] = (typeof(v) <: Union{AbstractVector, AbstractInterval}) ? v :
                     typeof(v) <: Integer ? (v:v) :
                     typeof(v) <: Number ? (v .. v) :
                     error("Bad index.")
    return OmniSpectrum(getindex(AxisArray(A); kwargs_new...))
end

function Base.setindex!(A::OmniSpectrum{US, UA} where {US, UA}, v; kwargs...)
    B = AxisArray(A)
    # size(v) == (1,1) && (v = v[1, 1]) # scalars
    setindex!(B, v; kwargs...)
    A = OmniSpectrum(B)
    return nothing
end

function Base.getindex(A::OmniSpectrum{US, UA} where {US, UA}, I::Union{Quantity, AbstractInterval})
    name = axistype(A.axis)
    kwargs = Dict()
    kwargs[name] = I
    return getindex(A; kwargs...)
end

function Base.setindex!(A::OmniSpectrum{US, UA} where {US, UA}, v, I::Union{Quantity, AbstractInterval})
    name = axistype(A.axis)
    kwargs = Dict()
    kwargs[name] = I
    B = AxisArray(A)
    setindex!(B, v; kwargs...)
    A = OmniSpectrum(B)
    return nothing
end

# support functions
axistype() = Dict(
    :direction => ((:direction, :direction), 𝐀),
    :frequency => ((:temporal, :frequency), 𝐓^-1),
    :angular_frequency => ((:temporal, :angular_frequency), 𝐀*𝐓^-1),
    :period => ((:temporal, :period), 𝐓),
    :angular_period => ((:temporal, :angular_period), 𝐓*𝐀^-1),
    :wavenumber => ((:spatial, :frequency), 𝐋^-1),
    :angular_wavenumber => ((:spatial, :angular_frequency), 𝐀*𝐋^-1),
    :wavelength => ((:spatial, :period), 𝐋),
    :angular_wavelength => ((:spatial, :angular_period), 𝐋*𝐀^-1)
)

axistype(ax::Symbol) = axistype()[ax]
axistype(dim::Dimensions)::Symbol = Dict(value[2] => key for (key, value) in axistype())[dim]
axistype(x::AbstractArray{<:Quantity}) = axistype(dimension(eltype(x)))
axistype(x::Quantity) = axistype(dimension(x))

function axistype(domain::Symbol, type::Symbol)::Symbol
    return Dict(value[1] => key for (key, value) in axistype())[(domain, type)]
end

function _check_typeconsistency(array::AbstractArray)::Nothing
    consistent = all(y->typeof(y)==eltype(array), array)
    !consistent && error("All elements of array must be of same type.")
    return nothing
end
