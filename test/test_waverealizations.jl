using Test
using AxisArrays
using Random
using WaveSpectra
using WaveRealizations

let
    amplitudes = ComplexAmplitudes(ones(2, 1), [1, 2] .* Hz, [90] .* °)
    @test WaveSpectra.coordinates(amplitudes) == :polar
    @test WaveSpectra.ispolar(amplitudes)
    @test !WaveSpectra.iscartesian(amplitudes)
end

let
    kx = [-1, 1] .* (rad / m)
    ky = [-2, 0, 2] .* (rad / m)
    amplitudes = ComplexAmplitudes(ones(2, 3), kx, ky)
    @test WaveSpectra.iscartesian(amplitudes)
end