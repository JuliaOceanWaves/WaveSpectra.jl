# TODO: Display of Spectrum and OmnidirectionalSpectrum in Pluto... combination of DataFrames and AxisArrays
# TODO: Option to display as images in html?

function Base.show(io::IO, ::MIME"text/html", s::Spectrum) # make rounding more flexible
    spectrum_size = join(size(s.data), "x")
    spectrum_coord = s.coordinates
    spectrum_unit = Spectra.unit(s, :spectrum)
    axis1_name, axis2_name = Spectra.axesnames(s)
    axis1_unit = Spectra.unit(s, :axis1)
    axis2_unit = Spectra.unit(s, :axis2)
    println(io, L"%$spectrum_size %$spectrum_coord spectrum (%$spectrum_unit) with axis1 = %$axis1_name (%$axis1_unit), axis2 = %$axis2_name (%$axis2_unit)")
    println(io)
    data = round.(Array(ustrip(s.data)), digits=2)
    Base.print_matrix(stdout, data)
end

function Base.show(io::IO, ::MIME"text/plain", s::Spectrum)
    spectrum_size = join(size(s.data), "x")
    spectrum_coord = s.coordinates
    spectrum_unit = Spectra.unit(s, :spectrum)
    axis1_name, axis2_name = Spectra.axesnames(s)
    axis1_unit = Spectra.unit(s, :axis1)
    axis2_unit = Spectra.unit(s, :axis2)
    println(io, "$spectrum_size $spectrum_coord spectrum ($spectrum_unit) with axis1 = $axis1_name ($axis1_unit), axis2 = $axis2_name ($axis2_unit)")
    println(io)
    data = round.(Array(ustrip(s.data)), digits=2)
    Base.print_matrix(stdout, data)
end

function Base.show(io::IO, s::Spectrum)
    if get(io, :compact, true)::Bool
        println(io, L"%$spectrum_size %$spectrum_coord spectrum (%$spectrum_unit) with axis1 = %$axis1_name (%$axis1_unit), axis2 = %$axis2_name (%$axis2_unit)")
    end
end