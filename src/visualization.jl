# TODO: investigate contour plot for polar spectrum

@recipe function f(x::Spectrum)
    if x.coordinates == :polar
        # data
        θ = ustrip.(uconvert.(rad, x.axis2))
        r = ustrip.(x.axis1)
        z = ustrip.(x.data)

        # plot info
        freq_name, angle_name = axesnames(x)
        freq_unit_str = repr(unit(x, :axis1), context=:fancy_exponent => true)
        # angle_unit = unit(x, :axis2)
        # angle_unit_str = repr(angle_unit, context=:fancy_exponent => true)
        spectrum_unit_str = repr(unit(x), context=:fancy_exponent => true)
        # angles = (0:π/4:(2π-π/4))*rad
        # angles_labels = map(angles) do angle
        #     repr(uconvert(angle_unit, angle), context=:fancy_exponent => true)
        # end

        # default properties
        seriestype --> :heatmap
        projection --> :polar
        xaxis --> false
        yaxis --> false
        label --> ""
        color --> :deep
        colorbar_title --> "Spectral Density\n($spectrum_unit_str)"
        right_margin --> 10*plots_mm

        # annotations
        n = 8
        z_pos = 1 * exp.(im * 2π * (0:n-1) / n)
        angle_unit = unit(x, :axis2)

        angles = (1τ |> angle_unit) * ((0:n-1) / n)
        angle_labels = map(angles) do angle
            label_str = repr(uconvert(angle_unit, angle), context=:fancy_exponent => true)
            return text(label_str, 12, "Computer Modern")
        end
        annotation = [(real(z_pos[i]), imag(z_pos[i]), angle_labels[i]) for i in 1:n]
        # annotation --> (real.(z), imag.(z_pos), text.(angle_labels, 12, "Computer Modern"))
        # annotation --> (-1.2, 0.9, text("Angular: $angle_name ($angle_unit)\nRadial: $freq_name ($freq_unit)", 8))
        annotation --> annotation
        return θ, r, z
    elseif x.coordinates == :cartesian
        freq1_name, freq2_name = axesnames(x)
        freq1_unit = unit(x, :axis1)
        freq2_unit = unit(x, :axis2)
        spectrum_unit = unit(x)

        seriestype --> :contour # use contour as defaul
        xlabel --> "$freq1_name ($freq1_unit)"
        ylabel --> "$freq2_name ($freq2_unit)"
        color --> :viridis
        colorbar_title --> "Spectral Density ($spectrum_unit)"
        fill --> true
        linewidth --> 0
        levels --> 100
        title --> "Cartesian Spectrum"

        return ustrip.(x.data)
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
