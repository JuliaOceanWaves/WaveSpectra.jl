using AxisArrays
using DimensionfulAngles.DefaultSymbols
using Integrals: SampledIntegralProblem, solve
using Test
using WaveSpectra

# Polar spectrum with evenly spaced axes
f = (1:5) .* 0.1Hz
θ = (0:3) .* 90°
data_polar = reshape(collect(1:(length(f) * length(θ))), length(f), length(θ)) .*
             (m^2 / Hz / °)
S_polar = Spectrum(data_polar, f, θ)

# Integrating over direction returns the omnidirectional spectrum.
omni = integrate(S_polar, :axis2)
manual_omni_data = ((data_polar[:, 1] .+ data_polar[:, end]) .* (step(θ) / 2)) .+
                   vec(sum(data_polar[:, 2:(end - 1)]; dims = 2)) .* step(θ)
@test omni isa OmnidirectionalSpectrum
@test omni.data ≈ manual_omni_data
@test omni ≈ integrate(S_polar, :direction)
@test omni ≈ OmnidirectionalSpectrum(S_polar)

# Integrating over frequency
directional = @test_logs (:warn, r"function of direction") integrate(S_polar, :axis1)
manual_directional = ((data_polar[1, :] .+ data_polar[end, :]) .* (step(f) / 2)) .+
                     vec(sum(data_polar[2:(end - 1), :]; dims = 1)) .* step(f)
@test directional isa AxisArray
@test axisvalues(directional) == (θ,)
@test directional ≈ manual_directional
@test directional ≈
      [integrate(OmnidirectionalSpectrum(data_polar[:, j], f)) for j in eachindex(θ)]
@test directional ==
      @test_logs (:warn, r"function of direction") integrate(S_polar, :frequency)

# Whole-spectrum integration matches integrating the omnidirectional spectrum
@test integrate(S_polar) ≈ integrate(omni)
@test integrate(S_polar; include_zero = true) ≈ integrate(omni; include_zero = true)

# Integrating over subrange of frequency axis
omni_idx = integrate(S_polar, :axis2, (2, 3))
omni_qty = integrate(S_polar, :axis2, (θ[2], θ[3]))
@test omni_idx isa OmnidirectionalSpectrum
@test omni_idx ≈ omni_qty
@test omni_idx == OmnidirectionalSpectrum(omni_idx.data, f)
@test omni_qty == OmnidirectionalSpectrum(omni_qty.data, f)

# Integrating over subrange of direction axis
directional_idx = @test_logs (:warn, r"function of direction") match_mode=:any integrate(
    S_polar, :axis1, (2, 4))
directional_qty = @test_logs (:warn, r"function of direction") match_mode=:any integrate(
    S_polar, :axis1, (f[2], f[4]))
@test directional_idx isa AxisArray
@test axisvalues(directional_idx) == (θ,)
@test directional_idx ≈ directional_qty

# Double integration over subranges
omni_sub = integrate(S_polar, :axis2, (θ[2], θ[4]))
@test integrate(S_polar, (f[2], f[4]), (θ[2], θ[4])) ≈ integrate(omni_sub, (f[2], f[4]))

# Invalid integrations throw error
@test_throws ArgumentError integrate(S_polar, :bad_axis)
@test_throws ArgumentError integrate(S_polar, :axis1, (10Hz, 11Hz))

# Omnidirectional spectrum
omni_manual = OmnidirectionalSpectrum([1, 2, 3, 4, 5] .* (m^2 / Hz), f)
@test integrate(omni_manual, (2, 4)) ≈ integrate(omni_manual, (f[2], f[4]))
axis_zero = vcat(0.0Hz, f)
data_zero = vcat(0.0m^2 / Hz, omni_manual.data)
@test integrate(omni_manual; include_zero = true) ≈
      integrate(AxisArray(data_zero, Axis{:frequency}(axis_zero)))
directional_zero = @test_logs (:warn, r"function of direction") integrate(
    S_polar, :axis1; include_zero = true)
@test directional_zero isa AxisArray
@test directional_zero ≈ [integrate(
           AxisArray(vcat(0.0m^2 / Hz / °, data_polar[:, j]), Axis{:frequency}(axis_zero))
       ) for j in eachindex(θ)]

# Cartesian spectrum
kx = collect((-2:2) .* (0.2 * rad / m))
ky = collect((-1:1) .* (0.5 * rad / m))
data_cart = reshape(collect(1:(length(kx) * length(ky))), length(kx), length(ky)) .*
            (m^4 / rad^2)
S_cart = Spectrum(data_cart, kx, ky)

spec_x = @test_logs (:warn, r"marginal spectrum") match_mode=:any integrate(
    S_cart, :axis2)
spec_y = @test_logs (:warn, r"marginal spectrum") match_mode=:any integrate(
    S_cart, :axis1)
@test spec_x isa AxisArray
@test spec_y isa AxisArray
@test axisvalues(spec_x) == (kx,)
@test axisvalues(spec_y) == (ky,)
@test spec_x == @test_logs (:warn, r"marginal spectrum") match_mode=:any integrate(
    S_cart, :angular_wavenumber_2)
@test spec_y == @test_logs (:warn, r"marginal spectrum") match_mode=:any integrate(
    S_cart, :angular_wavenumber_1)
@test integrate(S_cart) ≈ integrate(spec_x)
spec_x_sub = @test_logs (:warn, r"marginal spectrum") match_mode=:any integrate(
    S_cart, :axis2, (ky[1], ky[2]))
@test integrate(S_cart, (kx[2], kx[4]), (ky[1], ky[2])) ≈
      integrate(spec_x_sub, (kx[2], kx[4]))

# Rectangular integration
rect_omni = integrate(omni_manual; method = RectangularRule())
rect_range = integrate(omni_manual, (2, 4); method = RectangularRule())
rect_solution = solve(SampledIntegralProblem(omni_manual.data, f), RectangularRule())
@test rect_omni ≈ sum(omni_manual.data) * step(f)
@test rect_range ≈ sum(omni_manual.data[2:4]) * step(f)
@test rect_solution.u ≈ rect_omni
@test integrate(S_polar, :axis2; method = RectangularRule()).data ≈
      vec(sum(data_polar; dims = 2) .* step(θ))
weights_range = WaveSpectra.find_weights(1:4, RectangularRule())
weights_vec = WaveSpectra.find_weights(collect(f), RectangularRule())
@test weights_range isa WaveSpectra.RectangularUniformWeights
@test weights_vec isa WaveSpectra.RectangularUniformWeights
@test weights_range[1] == 1
@test weights_vec[1] == step(f)
@test_throws ArgumentError WaveSpectra.find_weights([1, 2, 4], RectangularRule())
nonuniform_omni = OmnidirectionalSpectrum([1.0, 2.0, 3.0] .* (m^2 / Hz),
    [0.1, 0.2, 0.4] .* Hz)
@test_throws ArgumentError integrate(nonuniform_omni; method = RectangularRule())
