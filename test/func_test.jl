using DimensionfulAngles: radᵃ as rad
using Unitful: Hz, m
using WaveSpectra

ex_units = [upreferred.(r) for r in keys(WaveSpectra._grad)]

for test_convert in ex_units
    println(test_convert)
    aux_spec = OmnidirectionalSpectrum(x -> x, typeof(1.0test_convert[1]))
    aux2_spec = convert_frequency(aux_spec, typeof(1.0test_convert[2]))
    
end