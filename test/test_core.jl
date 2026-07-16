
using AxisArrays
using DimensionfulAngles.DefaultSymbols
using Test
using WaveSpectra

# Polar spectra with temporal frequencies
Nf = 20
Δf = 0.1Hz
f = (1:Nf) * Δf
Nθ = 36
Δθ = 360° / Nθ
θ = (0:(Nθ - 1)) * Δθ

data = reshape(collect(1:(Nf * Nθ)), Nf, Nθ) .* (m^2 / Hz / °)
S = Spectrum(reverse(reverse(data; dims = 1); dims = 2), reverse(f), reverse(θ))

# Basic construction & properties.
@test size(S) == (Nf, Nθ)
@test eltype(S) == eltype(data)
@test WaveSpectra.ispolar(S)
@test !WaveSpectra.iscartesian(S)
@test WaveSpectra.coordinates(S) == :polar
@test WaveSpectra.axestypes(S) == (:frequency, :direction)
@test WaveSpectra.axesnames(S) == (:frequency, :direction)
@test all(WaveSpectra.AxisArrays.axes(S) .≈ (f, θ))
@test WaveSpectra.axesinfo(:frequency)[1] == (:temporal, :linear, :frequency)
@test WaveSpectra.axestypes(:temporal, :linear, :frequency) == :frequency
@test WaveSpectra.axesinfo(S) == (
    WaveSpectra.axesinfo(:frequency), WaveSpectra.axesinfo(:direction)
)
@test WaveSpectra.istemporal(f)
@test !WaveSpectra.isspatial(f)
@test WaveSpectra.islinear(f)
@test !WaveSpectra.isangular(f)
@test WaveSpectra.isfrequency(f)
@test !WaveSpectra.isperiod(f)
@test WaveSpectra.isspectralvariable(f)
@test WaveSpectra.isdirection(θ)
@test !WaveSpectra.isspectralvariable(θ)
@test !WaveSpectra.istemporal(θ)
@test !WaveSpectra.isspatial(θ)
@test WaveSpectra.isevenlyspaced(S)
let spacing = WaveSpectra.evenspacing(S)
    @test spacing[1][1] ≈ first(f)
    @test spacing[1][2] ≈ Δf
    @test spacing[1][3] == length(f)
    @test spacing[2][1] ≈ first(θ)
    @test spacing[2][2] ≈ Δθ
    @test spacing[2][3] == length(θ)
end
@test_throws ArgumentError Spectrum(data, (0:(Nf - 1)) * Δf, θ)
@test_throws ArgumentError Spectrum(data, ((-1):(Nf - 2)) * Δf, θ)
@test_throws ArgumentError Spectrum(data, vcat(f[1:(end - 1)], Inf * Hz), θ)

# Indexing by integer, axis name, and physical intervals.
s_11 = S[1, 1]
@test s_11 == S.data[1, 1]
s_named = S[frequency = 3, direction = 4]
@test s_named isa Spectrum
@test size(s_named) == (1, 1)
@test axisvalues(s_named) == (f[3:3], θ[4:4])
s_range = S[f[3] .. f[5], θ[4] .. θ[6]]
@test s_range isa Spectrum
@test size(s_range) == (3, 3)
@test axisvalues(s_range) == (f[3:5], θ[4:6])
s_range_spec = Spectrum(WaveSpectra.AxisArray(s_range))
@test s_range_spec.data == s_range.data
@test s_range_spec.axis1 == s_range.axis1
@test s_range_spec.axis2 == s_range.axis2

# Copy/similar.
S_copy = copy(S)
@test S_copy isa Spectrum
@test S_copy !== S
@test S_copy == S
@test similar(S, Float64, size(S)) isa Spectrum
@test similar(S, Float64, (2, 3)) isa Matrix{Float64}

# Broadcasting.
S_sum = S + S
@test S_sum isa Spectrum
@test S_sum == Spectrum(2 * S.data, S.axis1, S.axis2)
@test_throws DimensionMismatch S .+ Spectrum(data, f, θ .+ 0.1Δθ)
S_lt = S .< (10 * m^2 / Hz / °)
@test S_lt isa AbstractMatrix{Bool}
@test S_lt == (S.data .< (10 * m^2 / Hz / °))
@test (S .== S.data[1, 1]) == (S.data .== S.data[1, 1])
@test_throws DimensionMismatch S .< Spectrum(data, f, θ .+ 0.1Δθ)

# Conversion to/from `AxisArray`.
S_axis = WaveSpectra.AxisArray(S)
@test Spectrum(S_axis) == S

# Omnidirectional spectrum from a polar spectrum
omni = OmnidirectionalSpectrum(S)
@test omni isa OmnidirectionalSpectrum
@test length(omni) == Nf
@test WaveSpectra.axestypes(omni) == :frequency
@test WaveSpectra.axesnames(omni) == :frequency
@test all(WaveSpectra.AxisArrays.axes(omni) .≈ (f,))
@test WaveSpectra.isevenlyspaced(omni)
let spacing = WaveSpectra.evenspacing(omni)
    @test spacing[1] ≈ first(f)
    @test spacing[2] ≈ Δf
    @test spacing[3] == length(f)
end
@test_throws ArgumentError OmnidirectionalSpectrum(omni.data, (0:(Nf - 1)) * Δf)
@test_throws ArgumentError OmnidirectionalSpectrum(omni.data, ((-1):(Nf - 2)) * Δf)
@test_throws ArgumentError OmnidirectionalSpectrum(
    omni.data, vcat(f[1:(end - 1)], Inf * Hz)
)
@test_throws ArgumentError OmnidirectionalSpectrum(
    fill(1.0m^2 / s, 3), [0.0, 1.0, 2.0] .* s
)
omni_2 = omni[frequency = 2]
@test omni_2 isa AxisArray
@test size(omni_2) == (1,)
@test axisvalues(omni_2) == (f[2:2],)
@test OmnidirectionalSpectrum(omni_2) == OmnidirectionalSpectrum(omni.data[2:2], f[2:2])
@test OmnidirectionalSpectrum(WaveSpectra.AxisArray(omni)) == omni
@test omni[1] == omni.data[1]
@test omni[1] < (2 * omni.data[1])
omni_le = omni .<= (sum(omni.data) / length(omni.data))
@test omni_le isa AbstractVector{Bool}
@test omni_le == (omni.data .<= (sum(omni.data) / length(omni.data)))
@test (omni .== omni.data[1]) == (omni.data .== omni.data[1])
@test_throws DimensionMismatch omni .>= OmnidirectionalSpectrum(omni.data, f .+ 0.1Hz)

# Show methods.
@test occursin("Spectrum", sprint(show, S))
@test occursin("OmnidirectionalSpectrum", sprint(show, omni))
@test occursin("with polar coordinates and axes", sprint(show, MIME"text/plain"(), S))
@test occursin("with axis", sprint(show, MIME"text/plain"(), omni))

# Cartesian spectra with spatial frequencies:
# repeat some of the same tests above for this case.
Nkx = 8
Nky = 6
Δk = 0.25 * rad / m
kx = collect(((-(Nkx ÷ 2)):(Nkx ÷ 2)) * Δk)
ky = collect(((-(Nky ÷ 2)):(Nky ÷ 2)) * Δk)
data_k = reshape(collect(1:(length(kx) * length(ky))), length(kx), length(ky)) .*
         (m^4 / rad^2)
Sxy = Spectrum(data_k, reverse(kx), reverse(ky))
Sxy_inf = Spectrum(data_k, vcat(kx[1:(end - 1)], Inf * rad / m), ky)

@test size(Sxy) == (length(kx), length(ky))
@test isinf(Sxy_inf.axis1[end])
@test WaveSpectra.iscartesian(Sxy)
@test !WaveSpectra.ispolar(Sxy)
@test WaveSpectra.coordinates(Sxy) == :cartesian
@test WaveSpectra.axestypes(Sxy) == (:angular_wavenumber, :angular_wavenumber)
@test WaveSpectra.axesnames(Sxy) == (:angular_wavenumber_1, :angular_wavenumber_2)
@test all(WaveSpectra.AxisArrays.axes(Sxy) .≈ (kx, ky))
@test WaveSpectra.isspatial(kx)
@test !WaveSpectra.istemporal(kx)
@test WaveSpectra.isangular(kx)
@test !WaveSpectra.islinear(kx)
@test WaveSpectra.isfrequency(kx)
@test !WaveSpectra.isperiod(kx)
@test WaveSpectra.isspectralvariable(kx)
@test WaveSpectra.isevenlyspaced(Sxy)
sxy_ordered = Spectrum(data_k, kx, ky)
@test sxy_ordered.axis1 == kx
@test sxy_ordered.axis2 == ky
@test sxy_ordered.data == data_k
sxy_named = Sxy[angular_wavenumber_1 = 2, angular_wavenumber_2 = 3]
@test sxy_named isa Spectrum
@test size(sxy_named) == (1, 1)
@test axisvalues(sxy_named) == (kx[2:2], ky[3:3])
sxy_point = Sxy[kx[3], ky[4]]
@test sxy_point isa Spectrum
@test size(sxy_point) == (1, 1)
@test axisvalues(sxy_point) == (kx[3:3], ky[4:4])
@test Sxy[1, 1] == Sxy.data[1, 1]
sxy_range = Sxy[kx[3] .. kx[5], ky[2] .. ky[4]]
@test sxy_range isa Spectrum
@test size(sxy_range) == (3, 3)
@test axisvalues(sxy_range) == (kx[3:5], ky[2:4])
@test Spectrum(WaveSpectra.AxisArray(Sxy)) == Sxy
Sxy_copy = copy(Sxy)
@test Sxy_copy isa Spectrum
@test Sxy_copy !== Sxy
@test Sxy_copy == Sxy
@test similar(Sxy, Float64, size(Sxy)) isa Spectrum
@test similar(Sxy, Float64, (2, 3)) isa Matrix{Float64}

Sxy_sum = Sxy .+ Sxy
@test Sxy_sum isa Spectrum
@test Sxy_sum == Spectrum(2 .* Sxy.data, Sxy.axis1, Sxy.axis2)
@test_throws DimensionMismatch Sxy .+ Spectrum(data_k, kx, ky .+ 0.1Δk)
Sxy_gt = Sxy .> (10 * m^4 / rad^2)
@test Sxy_gt isa AbstractMatrix{Bool}
@test Sxy_gt == (Sxy.data .> (10 * m^4 / rad^2))

# axis types functions
T = (1:10) .* s
@test WaveSpectra.istemporal(T)
@test !WaveSpectra.isspatial(T)
@test WaveSpectra.islinear(T)
@test !WaveSpectra.isangular(T)
@test WaveSpectra.isperiod(T)
@test !WaveSpectra.isfrequency(T)
