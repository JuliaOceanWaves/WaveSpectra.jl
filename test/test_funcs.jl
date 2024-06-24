using Unitful: upreferred, uconvert, unit, dimension
using Unitful: 𝐓, Hz, s, 𝐋, m
using WaveSpectra
using Test

@testset "Spectra Functions" begin
    f = range(1.0Hz, 5.0Hz, 5)
    cont_spectrum = OmnidirectionalSpectrum(x->x*((m^2) / (Hz^2)), typeof(f[1]));
    disc_spectrum = DiscreteOmnidirectionalSpectrum(cont_spectrum, f);
    function test_spectrum(spectrum, args...)
        @test unit(spectrum) == m^2/Hz
        @test dimension(spectrum) == 𝐋^2*𝐓
        @test frequency_unit(spectrum) == Hz
        @test frequency_dimension(spectrum) == 𝐓^-1
        @test quantity(spectrum) == (𝐋^2, m^2)
        @test spectral_moment(spectrum, 0, args...) ≈ (12.0m^2)
        @test energy_period(spectrum, args...) ≈ (1s/3)
        @test significant_waveheight(spectrum, args...) ≈ (13.856406460551018m)
    end
    test_spectrum(cont_spectrum, 1.0Hz, 5.0Hz)
    test_spectrum(disc_spectrum)
end

ex_units = [upreferred.(r) for r in keys(WaveSpectra._grad_1)]
@testset "Conversion Equivalences Continous" begin
    for (old_TF, new_TF) in ex_units
        aux_spec = OmnidirectionalSpectrum(x -> x, typeof(1.0(old_TF)))
        new_spec = convert_frequency(aux_spec, 1.0(new_TF));

        int_begin, int_end = 1.0*(old_TF), 10.0*(old_TF)
        aux_begin = uconvert(new_TF, int_begin, WaveSpectra.deepwater_gradient)
        aux_end = uconvert(new_TF, int_end, WaveSpectra.deepwater_gradient)

        old_Hs = spectral_moment(aux_spec, 0, int_begin, int_end)
        new_Hs = spectral_moment(new_spec, 0, aux_begin, aux_end)
        @test old_Hs ≈ new_Hs
    end
end

@testset "Conversion Equivalences Discrete" begin
    for (old_TF, new_TF) in ex_units

        v=f=range(1.0, 5.0, 5)

        aux_spec = DiscreteOmnidirectionalSpectrum(v.*(old_TF), f.*(new_TF))
        new_spec = convert_frequency(aux_spec, 1.0(new_TF));

        old_Hs = spectral_moment(aux_spec, 0)
        new_Hs = spectral_moment(new_spec, 0)

        @test old_Hs ≈ new_Hs
    end
end
