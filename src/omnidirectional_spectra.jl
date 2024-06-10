

include("omnidirectional_spectra/continuous.jl")
include("omnidirectional_spectra/discrete.jl")


# # Plots recipes
# function _labels(spectrum::OmnidirectionalSpectrum)
#     x_label = "frequency"
#     y_label = isdensity(spectrum) ? "spectral density" : "discrete (integral) spectrum"
#     return (x_label, y_label)
# end

# @recipe function f(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}, args...) where {Ts, Tf, D}
#     _xlabel, _ylabel = _labels(spectrum)
#     xlabel --> _xlabel
#     ylabel --> _ylabel
#     if isdiscrete(spectrum)
#         marker := :auto
#         (spectrum.discrete.frequency, spectrum.discrete.value, args...)
#     else
#         isempty(args) && (args=(0*upreferred(Tf()), 10*upreferred(Tf())))
#         (spectrum.func, args...)
#     end
#     return nothing
# end
