
# Deep water waves
# module DeepWater
#     function steepness(spectrum::OmnidirectionalSpectrum{TS,TF},
#             f_begin::Union{Quantity,Nothing}=nothing, f_end::Union{Quantity,Nothing}=nothing;
#             alg::AbstractIntegralAlgorithm=QuadGKJL(), kwargs...) where {TS,TF}
#         hs = significant_waveheight(spectrum, f_begin, f_end; alg, kwargs...)
#         te = energy_period(spectrum, f_begin, f_end; alg, kwargs...)
#         return 2π * hs / (g * te ^ 2)
#     end
# end



# Parametric spectra
# te_to_tp(Te::Time, γ::Number) = Te / (0.8255 + 0.03852*γ - 0.005537*γ^2 + 0.0003154*γ^3)
# te_to_tp(Te::Time) = Te / 0.858

function pm_spectrum(significant_waveheight::Length, energy_period::Time)
    fp = uconvert(Hz, 0.858/energy_period)
    hs = significant_waveheight

    function spectrum(f::Frequency)
        if f==0Hz
            value = 0m^2/Hz
        else
            value = (hs^2 / 4) * (1.057 * fp)^4 * (f)^(-5) * exp((-5 / 4) * (fp / f)^4)
        end
        return value
    end
    return OmnidirectionalSpectrum(spectrum)
end

# pm_spectrum(Hs, Tp) = f -> pm_spectrum(Hs, Tp, f)

# function jonswap_spectrum(Hs, Tp, f; γ=nothing)
#     fp = 1 / Tp
#     σ = f <= fp ? 0.07 : 0.09
#     α = exp(-((f / fp - 1) / (√(2) * σ))^2)
#     isnothing(γ) && (
#         γ = (Tp / √Hs ≤ 3.6) ? 5 : (
#             (Tp / √Hs > 5) ? 1 : (
#                 exp(5.75 - 1.15(Tp / √Hs))
#             )))
#     Cws = 1 - 0.287 * log(γ)
#     spectrum = f == 0 ? 0 : Cws * pm_spectrum(Hs, Tp, f) * (γ^α)
#     return spectrum
# end

# jonswap_spectrum(Hs, Tp; γ=nothing) = f -> jonswap_spectrum(Hs, Tp, f; γ)
