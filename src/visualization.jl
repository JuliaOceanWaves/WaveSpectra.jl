# TODO: investigate contour plot for polar spectrum

@recipe function f(s::Spectrum)
    if s.coordinates == :polar

        freq = ustrip.(s.axis1)
        angle = ustrip.(s.axis2)
        x = deg2rad.(angle)
        y = freq
        z = s.data

        freq_name, angle_name = axesnames(s)
        freq_unit = unit(s, :axis1)
        angle_unit = unit(s, :axis2)
        spectrum_unit = unit(s)

        seriestype --> :heatmap
        projection --> :polar
        label --> ""
        color --> :viridis
        colorbar_title --> "Spectral Density ($spectrum_unit)"
        title --> "Polar Spectrum"
        annotation --> (-1.2, 0.9, text("Angular: $angle_name ($angle_unit)\nRadial: $freq_name ($freq_unit)", 8))

        return x, y, z
    elseif s.coordinates == :cartesian
        freq1_name, freq2_name = axesnames(s)
        freq1_unit = unit(s, :axis1)
        freq2_unit = unit(s, :axis2)
        spectrum_unit = unit(s)

        seriestype --> :contour # use contour as defaul
        xlabel --> "$freq1_name ($freq1_unit)"
        ylabel --> "$freq2_name ($freq2_unit)"
        color --> :viridis
        colorbar_title --> "Spectral Density ($spectrum_unit)"
        fill --> true
        linewidth --> 0
        levels --> 100
        title --> "Cartesian Spectrum"

        return ustrip.(s.data)
    end
end


@recipe function f(omni::OmnidirectionalSpectrum)

    x = omni.axis
    y = omni

    freq_name = axesnames(omni)

    xlabel --> "$freq_name"
    ylabel --> "Spectral density"
    label --> ""
    title --> "Omnidirectional Spectrum"

    return x, y
end