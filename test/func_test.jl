using DimensionfulAngles: radᵃ as rad
using Unitful: upreferred
using WaveSpectra
using Test

# TODO CHECK Convention for Unit Test Structs
# ex_units = [upreferred.(r) for r in keys(WaveSpectra._grad)]
# @testset "Conversions" begin 
#     for (old_TF, new_TF) in ex_units
#         @eval aux_spec = OmnidirectionalSpectrum(x -> x, typeof(1.0($old_TF)))
#         @test !isnothing(@eval convert_frequency(aux_spec, typeof(1.0($new_TF))))
#     end
# end