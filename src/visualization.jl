# TODO: add contour/surface plot using Makie and better 

@recipe function f(s::Spectrum)
    if s.coordinates == :polar
        if Spectra.isdirection(s.axis1)
            angle = ustrip.(s.axis1)
            freq = ustrip.(s.axis2)
            angle_name, freq_name = Spectra.axesnames(s)
            angle_unit = Spectra.unit(s, :axis1)
            freq_unit = Spectra.unit(s, :axis2)
        elseif Spectra.isdirection(s.axis2)
            freq = ustrip.(s.axis1)
            angle = ustrip.(s.axis2)
            freq_name, angle_name = Spectra.axesnames(s)
            freq_unit = Spectra.unit(s, :axis1)
            angle_unit = Spectra.unit(s, :axis2)
        end

        x = deg2rad.(repeat(angle, outer=length(freq)))
        y = repeat(freq, inner=length(angle))
        z = vec(ustrip.(s.data))
        spectrum_unit = Spectra.unit(s)

        seriestype := :scatter
        projection := :polar
        marker_z := z
        label := L""
        colorbar_title := L"Spectral Density (%$spectrum_unit)"
        title := L"Polar Spectrum"
        annotation := (-1.2, 0.9, text(L"Angular: %$angle_name (%$angle_unit)\nRadial: %$freq_name (%$freq_unit)", 8))

        x, y
    elseif s.coordinates == :cartesian # update for cartesian coordinates
        # For cartesian coordinates, just use a simple heatmap
        seriestype := :heatmap
        title := L"Cartesian Spectrum"

        Unitful.ustrip.(s.data)
    end
end


@recipe function f(omni::OmnidirectionalSpectrum)

    x = omni.axis
    y = omni

    freq_name = Spectra.axesnames(omni)

    xlabel := "$freq_name"
    ylabel := "Spectral density"
    title := "Omnidirectional Spectrum"

    x, y
end