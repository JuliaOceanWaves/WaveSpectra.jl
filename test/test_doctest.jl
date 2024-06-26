using Test, Documenter, WaveSpectra

ENV["UNITFUL_FANCY_EXPONENTS"] = true
ENV["JULIA_DEBUG"]=Documenter

doctest(WaveSpectra)

# Subexpression:
# │ 
# │ spectral_moment(s, 1, 1u"Hz", 5u"Hz")
# │ 
# │ Evaluated output:
# │ 
# │ 41.333333333333336 s⁻³
# │ 
# │ Expected output:
# │ 
# │ 41.333 s⁻³
# │ 
# │ 1 doctest filter was applied:
# │ 
# │   r"(\d*).(\d{3})\d+" => s"\1.\2"
# │ 
# │ Evaluated output (filtered):
# │ 
# │ 41.333333333333336 s⁻³
# │ 
# │ Expected output (filtered):
# │ 
# │ 41.333 s⁻³