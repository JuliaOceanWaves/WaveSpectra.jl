
function _check_omnidirectional_dimensions(TS, TF)
    if !(isa(1*upreferred(TF()), Frequency) || (TF == typeof(NoDims)))
        error("parameter `TF` must be a `Unitful.Frequency` or `Unitful.NoDims`")
    end
    if (TF == typeof(NoDims)) && !(TS == typeof(NoDims))
        error("if `TF` is dimensionless `TS` must be too")
    end
    return nothing
end

include("omnidirectional_spectra/continuous.jl")
