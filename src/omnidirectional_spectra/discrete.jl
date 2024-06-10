
# Struct
struct DiscreteOmnidirectionalSpectrum{S,F,D,N} <: AbstractArray{T, N} where {T}
    value::AbstractVecOrMat{<:Number}
    frequency::AbstractVector{<:Number}

    function DiscreteOmnidirectionalSpectrum(value, frequency; density=true)
        # parameters
        D = density
        N = ndims(value)
        TS = typeof(dimension(eltype(value)))
        TF = typeof(dimension(eltype(frequency)))
        # check parameters
        !isa(D, Bool) && error("parameter `D` must be a boolean")
        !(N∈[1,2]) && error("parameter `N` must be `1` or `2`")
        _check_omnidirectional_dimensions(TS, TF)
        # check arguments
        if (size(value)[1] ≠ length(frequency))
            error("`frequency` and `value` arrays must have same size")
        end
        # create parametric object
        return new{S,F,D,N}(value, frequency)
    end
end
