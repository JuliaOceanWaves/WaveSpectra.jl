# TODO: Display of Spectrum and OmnidirectionalSpectrum in Pluto... combination of DataFrames and AxisArrays
# TODO: Option to display as images in html?


# Compact show method
function Base.show(io::IO, x::Spectrum)
    shape = size(x.data);
    parameters = typeof(x).parameters
    us = unit(parameters[1])
    u1 = unit((parameters[2]).parameters[1])
    u2 = unit((parameters[3]).parameters[1])
    io_fancy = IOContext(io, :fancy_exponent => true)
    print(io_fancy, shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}")
end

# Detailed show method (text/plain)
function Base.show(io::IO, ::MIME"text/plain", x::Spectrum)
    shape = size(x.data)
    us, ui, u1, u2 = unit(x), unit(x, :integral), unit(x, :axis1), unit(x, :axis2)
    axt = axestypes(x)
    axt1 = titlecase(replace(String(axt[1]), "_" => " "))
    axt2 = titlecase(replace(String(axt[2]), "_" => " "))
    xc = String(x.coordinates)
    io_fancy = IOContext(io, :fancy_exponent => true)
    println(io_fancy, shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}",
        "\nSpectral density for Quantity (", ui,") with ", xc, " coordinates:",
        "\n  • Axis 1: ", axt1, " (", u1, ")\n  • Axis 2: ", axt2, " (", u2, ")",
        "\nand data:"
    )
    Base.print_matrix(io, ustrip.(x.data))
end

function Base.show(io::IO, ::MIME"text/html", x::Spectrum)
    shape = size(x.data)
    us = repr(unit(x), context=:fancy_exponent => true)
    ui = repr(unit(x, :integral), context=:fancy_exponent => true)
    u1 = repr(unit(x, :axis1), context=:fancy_exponent => true)
    u2 = repr(unit(x, :axis2), context=:fancy_exponent => true)
    axt = axestypes(x)
    axt1 = titlecase(replace(String(axt[1]), "_" => " "))
    axt2 = titlecase(replace(String(axt[2]), "_" => " "))
    xc = String(x.coordinates)
    println(io,
        "<p style='font-size:16px;'>",
        shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}",
        "<br>Spectral density for Quantity (", ui, ") with ", xc, " coordinates:",
        "<br>  • Axis 1: ", axt1, " (", u1, ")<br>  • Axis 2: ", axt2, " (", u2, ")",
        "<br>and data:</p>"
    )
    # title = ""
    # subtitle = ""
    pretty_table(
        io,
        ustrip.(x.data),
        backend=:html,
        column_labels = x.axis2,
        row_labels=x.axis1,
        maximum_number_of_rows=10,
        vertical_crop_mode=:middle,
        # title=title,
        # subtitle=subtitle,
        # title_alignment=:l,
        # subtitle_alignment=:l,
        style = HtmlTableStyle(;
            first_line_column_label = ["color" => "LightSeaGreen", "font-weight" => "bold"],
            row_label = ["color" => "Coral", "font-weight" => "bold"],
        )
    )
end

# function Base.show(io::IO, ::MIME"text/html", x::Spectrum)
#     shape = size(x.data)
#     us, ui, u1, u2 = unit(x), unit(x, :integral), unit(x, :axis1), unit(x, :axis2)
#     axt = axestypes(x)
#     axt1 = titlecase(replace(String(axt[1]), "_" => " "))
#     axt2 = titlecase(replace(String(axt[2]), "_" => " "))
#     xc = String(x.coordinates)
#     axis1 = x.axis1
#     axis2 = x.axis2
#     data = x.data
#     print(io, "<div class=\"wavespectra-spectrum\">",
#         "<style>",
#         ".wavespectra-spectrum{font-family:system-ui,Segoe UI,Helvetica,Arial,sans-serif;}",
#         ".wavespectra-spectrum table{border-collapse:collapse;font-size:0.9em;}",
#         ".wavespectra-spectrum th,.wavespectra-spectrum td{border:1px solid #ccc;padding:2px 6px;text-align:right;}",
#         ".wavespectra-spectrum th{background:#f5f5f5;font-weight:600;}",
#         ".wavespectra-spectrum .summary{margin-bottom:6px;}",
#         "</style>",
#         "<div class=\"summary\">",
#         shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}",
#         "<br>Spectral density for Quantity (", ui, ") with ", xc, " coordinates:",
#         "<br>Axis 1: ", axt1, " (", u1, "), Axis 2: ", axt2, " (", u2, ")",
#         "</div>",
#         "<table><thead><tr><th>", axt1, " \\ ", axt2, "</th>"
#     )
#     @inbounds for j in eachindex(axis2)
#         print(io, "<th>", ustrip(axis2[j]), "</th>")
#     end
#     print(io, "</tr></thead><tbody>")
#     @inbounds for i in eachindex(axis1)
#         print(io, "<tr><th>", ustrip(axis1[i]), "</th>")
#         @inbounds for j in eachindex(axis2)
#             print(io, "<td>", ustrip(data[i, j]), "</td>")
#         end
#         print(io, "</tr>")
#     end
#     print(io, "</tbody></table></div>")
# end
