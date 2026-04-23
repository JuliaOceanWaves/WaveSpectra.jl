module WaveSpectraMakieExt

import WaveSpectra
using Unitful: uconvert, unit, ustrip
using DimensionfulAngles: radᵃ as rad, turnᵃ as τ
using Makie: Axis as MAxis, Colorbar, Figure, PolarAxis, contourf!, lines!

function WaveSpectra.plot_spectrum!(ax, s::WaveSpectra.AbstractSpectrum; (func!) = contourf!,
        make_cyclic = true, ax_properties = Dict(), kwargs...)
    for (prop, val) in ax_properties
        setproperty!(ax, prop, val)
    end

    function setaxproperty!(prop, val)
        (prop ∉ keys(ax_properties)) && setproperty!(ax, prop, val)
    end

    if WaveSpectra.ispolar(s)
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

        θ = ustrip.(uconvert.(rad, s.axis2))
        r = collect(ustrip.(s.axis1))
        z = ustrip.(s.data)'
        if make_cyclic
            push!(θ, 2π)
            z = vcat(z, z[1:1, :])
        end

        plot_props = merge(Dict(:colormap => :deep), Dict(kwargs))
        p = func!(ax, θ, r, z; plot_props...)
    else
        x = ustrip.(s.axis1)
        y = ustrip.(s.axis2)
        z = ustrip.(s.data)

        freq1_name, freq2_name = WaveSpectra.axesnames(s)
        freq1_unit = repr(unit(s, :axis1), context = :fancy_exponent => true)
        freq2_unit = repr(unit(s, :axis2), context = :fancy_exponent => true)

        setaxproperty!(:xlabel, "$freq1_name ($freq1_unit)")
        setaxproperty!(:ylabel, "$freq2_name ($freq2_unit)")
        setaxproperty!(:xgridcolor, :seashell4)
        setaxproperty!(:ygridcolor, :seashell4)

        plot_props = merge(Dict(:colormap => :deep), Dict(kwargs))
        p = func!(ax, x, y, z; plot_props...)
    end
    return p
end

function WaveSpectra.plot_spectrum(s::WaveSpectra.AbstractSpectrum;
        (func!) = contourf!,
        make_cyclic = true,
        ax_properties = Dict(),
        colorbar_properties = Dict(),
        kwargs...)
    spectrum_unit = repr(unit(s), context = :fancy_exponent => true)
    colorbar_properties = merge(
        Dict(:label => "Spectral Density\n($spectrum_unit)"),
        colorbar_properties,
    )

    f = Figure()
    ax = WaveSpectra.ispolar(s) ? PolarAxis(f[1, 1]) : MAxis(f[1, 1])
    p = WaveSpectra.plot_spectrum!(ax, s; func!, make_cyclic, ax_properties, kwargs...)
    Colorbar(f[1, 2], p; colorbar_properties...)
    return f
end

function WaveSpectra.plot_spectrum!(ax, s::WaveSpectra.AbstractOmnidirectionalSpectrum;
        (func!) = lines!, ax_properties = Dict(), kwargs...)
    for (prop, val) in ax_properties
        setproperty!(ax, prop, val)
    end

    function setaxproperty!(prop, val)
        (prop ∉ keys(ax_properties)) && setproperty!(ax, prop, val)
    end

    x = ustrip.(s.axis)
    y = ustrip.(s.data)

    freq_name = WaveSpectra.axesnames(s)
    freq_unit = repr(unit(s, :axis), context = :fancy_exponent => true)
    spectrum_unit = repr(unit(s), context = :fancy_exponent => true)

    setaxproperty!(:xlabel, "$freq_name ($freq_unit)")
    setaxproperty!(:ylabel, "Spectral Density ($spectrum_unit)")
    setaxproperty!(:limits, ((0, nothing), (0, nothing)))

    plot_props = merge(Dict(), Dict(kwargs))
    return func!(ax, x, y; plot_props...)
end

function WaveSpectra.plot_spectrum(s::WaveSpectra.AbstractOmnidirectionalSpectrum;
        (func!) = lines!, ax_properties = Dict(), kwargs...)
    f = Figure()
    ax = MAxis(f[1, 1])
    WaveSpectra.plot_spectrum!(ax, s; func!, ax_properties, kwargs...)
    return f
end

end
