# plotting functions and recipes

# TODO: Make these package extensions for PLots and Makie?
# Plots.jl
@recipe function f(x::AbstractSpectrum)
    if x.coordinates == :polar
        # data
        θ = ustrip.(uconvert.(rad, x.axis2))
        r = ustrip.(x.axis1)
        z = ustrip.(x.data)
        spectrum_unit = repr(unit(x), context = :fancy_exponent => true)

        seriestype --> :heatmap
        projection --> :polar
        label --> ""
        color --> :deep
        colorbar_title --> "Spectral Density\n($spectrum_unit)"
        right_margin --> 10 * plots_mm
        return θ, r, z
    elseif x.coordinates == :cartesian
        freq1_name, freq2_name = axesnames(x)
        freq1_unit = repr(unit(x, :axis1), context = :fancy_exponent => true)
        freq2_unit = repr(unit(x, :axis2), context = :fancy_exponent => true)
        spectrum_unit = repr(unit(x), context = :fancy_exponent => true)

        seriestype --> :contour
        xlabel --> "$freq1_name ($freq1_unit)"
        ylabel --> "$freq2_name ($freq2_unit)"
        color --> :deep
        colorbar_title --> "Spectral Density\n($spectrum_unit)"
        fill --> true
        linewidth --> 0
        levels --> 100
        return ustrip.(x.data)
    end
end

@recipe function f(x::AbstractOmnidirectionalSpectrum)
    _x = ustrip.(x.axis)
    _y = ustrip.(x.data)
    freq_name = axesnames(x)
    freq_unit = repr(unit(x, :axis), context = :fancy_exponent => true)
    spectrum_unit = repr(unit(x), context = :fancy_exponent => true)

    xlabel --> "$freq_name ($freq_unit)"
    ylabel --> "Spectral Density ($spectrum_unit)"
    label --> ""
    return _x, _y
end

