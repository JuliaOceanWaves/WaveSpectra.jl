using Test
using Plots: Plots
using Makie: Axis as MAxis, Figure, PolarAxis
using WaveSpectra

# Polar spectrum
f = (1:3) * 0.1Hz
θ = (0:3) * 90°
polar = Spectrum(reshape(1:12, 3, 4) .* (m^2 / Hz / °), f, θ)
omni = OmnidirectionalSpectrum(polar)

# Cartesian spectrum
k = collect((-1:1) * (0.2 * rad / m))
cartesian = Spectrum(reshape(1:9, 3, 3) .* (m^4 / rad^2), k, k)
cartesian_period = uconvert(m, :axis1, uconvert(m, :axis2, cartesian))

# Test Plots
@test Plots.plot(polar) isa Plots.Plot
@test Plots.plot(polar; include_zero = true) isa Plots.Plot
@test Plots.plot(cartesian) isa Plots.Plot
@test Plots.plot(cartesian_period) isa Plots.Plot
@test Plots.plot(omni) isa Plots.Plot
@test Plots.plot(omni; include_zero = true) isa Plots.Plot

# Test Makie: new figures
@test plot_spectrum(polar) isa Figure
@test plot_spectrum(polar; include_zero = true) isa Figure
@test plot_spectrum(cartesian) isa Figure
@test plot_spectrum(cartesian_period) isa Figure
@test plot_spectrum(omni) isa Figure
@test plot_spectrum(omni; include_zero = true) isa Figure

# Test Makie: plotting into existing axes
fig = Figure()
polar_ax = PolarAxis(fig[1, 1])
polar_plot = plot_spectrum!(
    polar_ax,
    polar;
    make_cyclic = false,
    include_zero = true,
    ax_properties = Dict(:rgridcolor => :red)
)
@test polar_plot !== nothing
@test polar_ax.rgridcolor[] == :red
wrong_ax = MAxis(fig[1, 2])
@test_throws ArgumentError plot_spectrum!(wrong_ax, polar)

# Test Makie: plotting cartesian spectra
cartesian_ax = MAxis(fig[2, 1])
cartesian_plot = plot_spectrum!(
    cartesian_ax,
    cartesian;
    ax_properties = Dict(:xlabel => "custom x")
)
@test cartesian_ax.xlabel[] == "custom x"
@test occursin("angular_wavenumber_2", cartesian_ax.ylabel[])

# Test Makie: omnidirectional spectrum
omni_ax = MAxis(fig[2, 2])
omni_plot = plot_spectrum!(
    omni_ax,
    omni;
    include_zero = true,
    ax_properties = Dict(:ylabel => "custom y")
)
@test occursin("frequency", omni_ax.xlabel[])
@test omni_ax.ylabel[] == "custom y"
