module WaveSpectraMakieExt

using DimensionfulAngles: turnᵃ as τ
using Makie: Axis as MAxis, Colorbar, Figure, PolarAxis, contourf!, lines!
using Unitful: uconvert, unit, ustrip
using WaveSpectra: AbstractOmnidirectionalSpectrum, AbstractSpectrum, axesnames, ispolar,
                   rad,
                   _prepend_zero_axis, _prepend_zero_data

import WaveSpectra: plot_spectrum, plot_spectrum!

"""
    plot_spectrum!(ax, s::AbstractSpectrum; func! = contourf!, make_cyclic = true,
                   include_zero = false, ax_properties = Dict(), kwargs...)
    plot_spectrum!(ax, s::AbstractOmnidirectionalSpectrum;
                   (func!)=lines!, include_zero = false, ax_properties=Dict(), kwargs...)

Plot a spectrum `s` into an existing Makie axis `ax`.

For polar spectra, `ax` must be a `PolarAxis`.

See [`plot_spectrum`](@ref) for more information.
"""
function plot_spectrum!(ax, s::AbstractSpectrum; (func!) = contourf!,
        make_cyclic = true, include_zero = false, ax_properties = Dict(), kwargs...)
    for (prop, val) in ax_properties
        setproperty!(ax, prop, val)
    end

    function setaxproperty!(prop, val)
        (prop ∉ keys(ax_properties)) && setproperty!(ax, prop, val)
    end

    if ispolar(s)
        !(ax isa PolarAxis) && throw(ArgumentError(
            "Axis must be `PolarAxis` for a polar `Spectrum`")
        )
        runit = unit(s, :axis1)
        runit_str = repr(runit, context = :fancy_exponent => true)
        θunit = unit(s, :axis2)
        thetaticks_pos = collect(0:(1 / 8):(7 / 8)) * 2π
        thetaticks_label = (
            (θunit == rad) ? ["0, 2π", "¼π", "½π", "¾π", "π", "1¼π", "1½π", "1¾π"] :
            ((θunit == τ) ? ["0, τ", "⅛τ", "¼τ", "⅜τ", "½τ", "⅝τ", "¾τ", "⅞τ"] :
             ["0°", "45°", "90°", "135°", "180°", "225°", "270°", "315°"])
        )

        setaxproperty!(:thetaticks, (thetaticks_pos, thetaticks_label))
        setaxproperty!(:rtickformat, values -> ["$(v) $(runit_str)" for v in values])
        setaxproperty!(:rgridcolor, :seashell4)
        setaxproperty!(:thetagridcolor, :seashell4)
        setaxproperty!(:rtickangle, π / 2)
        setaxproperty!(:gridz, 1)

        axis1 = include_zero ? _prepend_zero_axis(s.axis1) : s.axis1
        data = include_zero ? _prepend_zero_data(s.data) : s.data
        θ = ustrip.(uconvert.(rad, s.axis2))
        r = collect(ustrip.(axis1))
        z = ustrip.(data)'
        if make_cyclic
            push!(θ, 2π)
            z = vcat(z, z[1:1, :])
        end

        defaults = Dict(
            :colormap => :deep,
        )
        plot_props = merge(defaults, Dict(kwargs))
        p = func!(ax, θ, r, z; plot_props...)
    else
        x = ustrip.(s.axis1)
        y = ustrip.(s.axis2)
        z = ustrip.(s.data)

        freq1_name, freq2_name = axesnames(s)
        freq1_unit = repr(unit(s, :axis1), context = :fancy_exponent => true)
        freq2_unit = repr(unit(s, :axis2), context = :fancy_exponent => true)

        setaxproperty!(:xlabel, "$freq1_name ($freq1_unit)")
        setaxproperty!(:ylabel, "$freq2_name ($freq2_unit)")
        setaxproperty!(:xgridcolor, :seashell4)
        setaxproperty!(:ygridcolor, :seashell4)

        defaults = Dict(
            :colormap => :deep,
        )
        plot_props = merge(defaults, Dict(kwargs))
        p = func!(ax, x, y, z; plot_props...)
    end
    return p
end

function plot_spectrum!(ax, s::AbstractOmnidirectionalSpectrum;
        (func!) = lines!, include_zero = false, ax_properties = Dict(), kwargs...
)
    for (prop, val) in ax_properties
        setproperty!(ax, prop, val)
    end

    function setaxproperty!(prop, val)
        (prop ∉ keys(ax_properties)) && setproperty!(ax, prop, val)
    end

    axis = include_zero ? _prepend_zero_axis(s.axis) : s.axis
    data = include_zero ? _prepend_zero_data(s.data) : s.data
    x = ustrip.(axis)
    y = ustrip.(data)

    freq_name = axesnames(s)
    freq_unit = repr(unit(s, :axis), context = :fancy_exponent => true)
    spectrum_unit = repr(unit(s), context = :fancy_exponent => true)

    setaxproperty!(:xlabel, "$freq_name ($freq_unit)")
    setaxproperty!(:ylabel, "Spectral Density ($spectrum_unit)")
    setaxproperty!(:limits, ((0, nothing), (0, nothing)))

    defaults = Dict()

    plot_props = merge(defaults, Dict(kwargs))
    func!(ax, x, y; plot_props...)
end

"""
    plot_spectrum(s::AbstractSpectrum; func! = contourf!, make_cyclic = true,
                  include_zero = false, ax_properties = Dict(),
                  colorbar_properties = Dict(), kwargs...)
    plot_spectrum(s::AbstractOmnidirectionalSpectrum; func! = lines!,
                  include_zero = false, ax_properties = Dict(), kwargs...)

Create a new Makie figure and plot spectrum `s`.

For `AbstractSpectrum`, this creates a plot with a colorbar.
For `AbstractOmnidirectionalSpectrum`, this creates a line plot.
"""
function plot_spectrum(s::AbstractSpectrum;
        (func!) = contourf!,
        make_cyclic = true,
        include_zero = false,
        ax_properties = Dict(),
        colorbar_properties = Dict(),
        kwargs...
)
    spectrum_unit = repr(unit(s), context = :fancy_exponent => true)
    defaults = Dict(
        :label => "Spectral Density\n($spectrum_unit)",
    )
    colorbar_properties = merge(defaults, colorbar_properties)

    f = Figure()
    ax = ispolar(s) ? PolarAxis(f[1, 1]) : MAxis(f[1, 1])
    p = plot_spectrum!(ax, s; func!, make_cyclic, include_zero, ax_properties, kwargs...)
    Colorbar(f[1, 2], p; colorbar_properties...)
    return f
end

function plot_spectrum(s::AbstractOmnidirectionalSpectrum;
        (func!) = lines!, include_zero = false, ax_properties = Dict(), kwargs...
)
    f = Figure()
    ax = MAxis(f[1, 1])
    plot_spectrum!(ax, s; func!, include_zero, ax_properties, kwargs...)
    return f
end

end
