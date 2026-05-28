using Test
using WaveSpectra: rad, m
using WaveSpectra.DispersionRelations

k = 2π * rad / m  # wavelength = 1 m
deep_depth = 1.0m
shallow_depth = 0.01m

deepwater = gravitywaves_deepwater()
shallowwater = gravitywaves_shallowwater(shallow_depth)
finite_deep = gravitywaves(deep_depth)
finite_shallow = gravitywaves(shallow_depth)

ω_deep = deepwater.dispersion(k)
ω_shallow = shallowwater.dispersion(k)
ω_finite_deep = finite_deep.dispersion(k)
ω_finite_shallow = finite_shallow.dispersion(k)

# test roundtrip
@test deepwater.dispersion_inverse(ω_deep)≈k rtol=1e-10
@test shallowwater.dispersion_inverse(ω_shallow)≈k rtol=1e-10
@test finite_deep.dispersion_inverse(ω_finite_deep)≈k rtol=1e-10
@test finite_shallow.dispersion_inverse(ω_finite_shallow)≈k rtol=1e-10

# test deep/shallow water approximations
@test ω_finite_deep≈ω_deep rtol=5e-3
@test finite_deep.gradient(k)≈deepwater.gradient(k) rtol=5e-3
@test ω_finite_shallow≈ω_shallow rtol=5e-3
@test finite_shallow.gradient(k)≈shallowwater.gradient(k) rtol=5e-3
