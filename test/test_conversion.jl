using Test
using DimensionfulAngles.DefaultSymbols
using Unitful
using WaveSpectra

# Build representative omnidirectional and polar directional spectra
dispersion = WaveSpectra.DispersionRelations.gravitywaves_deepwater()
f = range(0.01, 0.5, length = 512) .* Hz
omni = WaveSpectra.ParametricSpectra.spectrum_pierson_moskowitz(f, 2.5m, 8s)
spread = WaveSpectra.ParametricSpectra.spread_cartwright(
    range(0, 2π, length = 73) .* rad,
    omni.axis,
    π / 3 * rad,
    π / 8 * rad
)
spec = Spectrum(omni, spread)

axis_units = (
    Hz,
    rad / s,
    s,
    s / rad,
    m^-1,
    rad / m,
    m,
    m / rad
)

omni_variance = integrate(omni)
spec_variance = integrate(spec)

# Convert to all 8 possible frequency types and check total variance remains constant
for u in axis_units
    omni_u = uconvert(u, :axis, omni, dispersion)
    spec_u = uconvert(u, :axis1, spec, dispersion)

    @test integrate(omni_u)≈omni_variance rtol=5e-2
    @test integrate(spec_u)≈spec_variance rtol=5e-2
    @test unit(omni_u, :axis) == u
    @test unit(spec_u, :axis1) == u
end

# Convert the directional axis
spec_deg = uconvert(°, :direction, spec)
@test integrate(spec_deg) ≈ spec_variance
@test unit(spec_deg, :axis2) == °
@test spec_deg.axis2 ≈ uconvert.(°, spec.axis2)

# Convert units of the integral quantity
@test unit(uconvert(unit(omni, :integral), omni), :integral) == unit(omni, :integral)
@test unit(uconvert(unit(spec, :integral), spec), :integral) == unit(spec, :integral)
@test unit(uconvert(unit(spec, :axis2), :axis2, spec), :axis2) == unit(spec, :axis2)
@test integrate(uconvert(cm^2, omni)) ≈ omni_variance
@test integrate(uconvert(cm^2, spec)) ≈ spec_variance

# Test errors thrown
@test_throws ArgumentError uconvert(unit(spec), :spectrum, spec)
@test_throws ArgumentError uconvert(unit(omni), :spectrum, omni)
@test_throws ArgumentError uconvert(Hz, :not_an_axis, omni)
@test_throws ArgumentError uconvert(Hz, :not_an_axis, spec)

# Cartesian spectrum
Ncart = 256
cartesian_axis = collect(((-(Ncart ÷ 2)):(Ncart ÷ 2)) .* (0.3 / (Ncart ÷ 2) * rad / m))
profile = (x -> sin(x) / 2 + 0.5).(range(3π / 2, 7π / 2, length = length(cartesian_axis)))
cartesian_data = profile .* profile'
cartesian_spec = Spectrum(cartesian_data, reverse(cartesian_axis), cartesian_axis)
cartesian_variance = integrate(cartesian_spec)

for u in (m^-1, rad / m, m, m / rad)
    cartesian_u = uconvert(u, :axis1, uconvert(u, :axis2, cartesian_spec))

    if WaveSpectra.isperiod(u)
        @test isinf(cartesian_u.axis1[end])
        @test isinf(cartesian_u.axis2[end])
    else
        @test integrate(cartesian_u)≈cartesian_variance rtol=5e-2
    end
    @test unit(cartesian_u, :axis1) == u
    @test unit(cartesian_u, :axis2) == u
end

# Polar/cartesian point-cloud conversions
polar_points = polar_to_cartesian(spec)
@test polar_points isa WaveSpectra.AxisArray
@test size(polar_points) == (length(spec.axis1) * length(spec.axis2), 3)
@test WaveSpectra.axisvalues(polar_points)[2] == [:frequency_x, :frequency_y, :spectrum]
@test polar_points[1, 1] == spec.axis1[1] * cos(spec.axis2[1])
@test polar_points[1, 2] == spec.axis1[1] * sin(spec.axis2[1])
@test polar_points[1, 3] == spec.data[1, 1] * WaveSpectra.θ₀ / spec.axis1[1]
@test_throws ArgumentError Spectrum(
    [1.0 1.0; 3.0 4.0] .* (m^2 / Hz / rad),
    [0.0, 0.1] .* Hz,
    [0.0, π / 2] .* rad
)

cartesian_points = cartesian_to_polar(cartesian_spec)
@test cartesian_points isa WaveSpectra.AxisArray
@test size(cartesian_points) == (length(cartesian_spec.axis1) * length(cartesian_spec.axis2), 3)
@test WaveSpectra.axisvalues(cartesian_points)[2] == [:angular_wavenumber, :direction, :spectrum]
@test cartesian_points[1, 1] == sqrt(cartesian_spec.axis1[1]^2 + cartesian_spec.axis2[1]^2)
@test cartesian_points[1, 2] == atan(°, cartesian_spec.axis2[1], cartesian_spec.axis1[1])
@test cartesian_points[1, 3] ==
      cartesian_spec.data[1, 1] * cartesian_points[1, 1] / WaveSpectra.θ₀
