using Unitful: upreferred, uconvert, unit, dimension, ustrip
using Unitful: 𝐓, Hz, s, 𝐋, m
using WaveSpectra
using Test

function box(x)
    ax = ustrip(x)
    if (ax == 0.5 || ax == 1.5)
        return 0.5
    elseif (0.5 < ax < 1.5)
        return 1.0
    else
        return 0.0
    end
end
f = range(0.1, 2.1, 51)

@testset "Spectra Functions" begin
    cont_spectrum = OmnidirectionalSpectrum(x -> box(ustrip(x)).*m^2/Hz, typeof(1.0Hz));
    disc_spectrum = DiscreteOmnidirectionalSpectrum(cont_spectrum, f.*Hz);
    function test_spectrum(spectrum; atol=1e-2)
        @test unit(spectrum) == m^2/Hz
        @test dimension(spectrum) == 𝐋^2*𝐓
        @test frequency_unit(spectrum) == Hz
        @test frequency_dimension(spectrum) == 𝐓^-1
        @test quantity(spectrum) == (𝐋^2, m^2)
        @test spectral_moment(spectrum, 0) ≈ (1.0m^2)  atol=atol*m^2
        @test energy_period(spectrum) ≈ (1.09s)         atol=atol*s
        @test significant_waveheight(spectrum) ≈ (4m)   atol=atol*m
    end
    test_spectrum(cont_spectrum)
    test_spectrum(disc_spectrum)
end

ex_units = [upreferred.(r) for r in keys(WaveSpectra._grad_1)]
@testset "Conversion Equivalences Continous" begin
    println("Discrete")
    for (old_TF, new_TF) in ex_units
        
        aux_spec = OmnidirectionalSpectrum(x -> box(ustrip(x)).*(m^2/(old_TF)), typeof(1.0(old_TF)))
        new_spec = convert_frequency(aux_spec, 1.0(new_TF));

        # @test_nowarn significant_waveheight(aux_spec)
        # @test_nowarn significant_waveheight(new_spec)
        old_int = spectral_moment(aux_spec, 0)
        new_int = spectral_moment(new_spec, 0)
        @test old_int ≈ new_int atol=1e-2m^2
        @test new_int ≈ 1.0m^2 atol=1e-2m^2
        # println(old_TF, "->", new_TF)
    end
end

@testset "Conversion Equivalences Discrete" begin
    for (old_TF, new_TF) in ex_units
        aux_spec = DiscreteOmnidirectionalSpectrum(x -> box(ustrip(x)).*(m^2/(old_TF)), f.*(old_TF))
        new_spec = convert_frequency(aux_spec, 1.0(new_TF));

        # old_int = significant_waveheight(aux_spec)
        # new_int = significant_waveheight(new_spec)
        # @test old_int ≈ new_int atol=1e-2m
        old_int = spectral_moment(aux_spec, 0)
        new_int = spectral_moment(new_spec, 0)
        @test old_int ≈ new_int atol=1e-2m^2
        @test new_int ≈ 1.0m^2 atol=1e-2m^2
        # println(old_TF, "->", new_TF)
    end
end
