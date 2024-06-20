using DimensionfulAngles: radᵃ as rad
using Unitful: upreferred, uconvert
using WaveSpectra
using Test

ex_units = [upreferred.(r) for r in keys(WaveSpectra._grad_1)]
@testset "Conversion Equivalences" begin 
    for (old_TF, new_TF) in ex_units
        aux_spec = OmnidirectionalSpectrum(x -> x, typeof(1.0(old_TF)))
        # @test_nowarn convert_frequency(aux_spec, 1.0(new_TF));
        new_spec = convert_frequency(aux_spec, 1.0(new_TF));

        int_begin, int_end = 1.0*(old_TF), 10.0*(old_TF)
        aux_begin = uconvert(new_TF, int_begin, WaveSpectra.deepwater_gradient)
        aux_end = uconvert(new_TF, int_end, WaveSpectra.deepwater_gradient)

        old_Hs = spectral_moment(aux_spec, 0, int_begin, int_end)
        new_Hs = spectral_moment(new_spec, 0, aux_begin, aux_end)
        # println("$(dimension(old_TF)) --> $(dimension(new_TF))")
        @test old_Hs ≈ new_Hs
    end
end