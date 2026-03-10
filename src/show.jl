
# Compact show method
function Base.show(io::IO, x::Spectrum)
    shape = size(x.data);
    us, u1, u2 = unit(x), unit(x, :axis1), unit(x, :axis2)
    (us isa FreeUnits{(), NoDims, nothing}) && (us = 1)
    io_fancy = IOContext(io, :fancy_exponent => true)
    print(io_fancy, shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}")
end

# Detailed show method (text/plain)
function Base.show(io::IO, ::MIME"text/plain", x::Spectrum)
    shape = size(x.data)
    us, ui, u1, u2 = unit(x), unit(x, :integral), unit(x, :axis1), unit(x, :axis2)
    (us isa FreeUnits{(),NoDims,nothing}) && (us = 1)
    axt = axestypes(x)
    axt1 = titlecase(replace(String(axt[1]), "_" => " "))
    axt2 = titlecase(replace(String(axt[2]), "_" => " "))
    xc = String(x.coordinates)
    io_fancy = IOContext(io, :fancy_exponent => true)
    println(io_fancy,
        shape[1], "×", shape[2], " Spectrum{", us, "}{", u1, "}{", u2, "}",
        "\nSpectral density for Quantity (", ui,") with ", xc, " coordinates:",
        "\n  • Axis 1: ", axt1, " (", u1, ")",
        "\n  • Axis 2: ", axt2, " (", u2, ")",
        "\nand data(", us, "):"
    )
    Base.print_matrix(io, ustrip.(x.data))
end

# HTML show method (Pluto, Jupyter)
function Base.show(io::IO, m::MIME"text/html", x::Spectrum)
    s = size(x)
    us = unit(x)
    (us isa FreeUnits{(),NoDims,nothing}) && (us = 1)
    us = repr(us, context=:fancy_exponent => true)
    ui = repr(unit(x, :integral), context=:fancy_exponent => true)
    u1 = repr(unit(x, :axis1), context=:fancy_exponent => true)
    u2 = repr(unit(x, :axis2), context=:fancy_exponent => true)
    axt = axestypes(x)
    axt1 = titlecase(replace(String(axt[1]), "_" => " "))
    axt2 = titlecase(replace(String(axt[2]), "_" => " "))
    xc = String(x.coordinates)
    println(io,
        "<p><strong>", s[1], "×", s[2], " Spectrum{", us, "}{", u1, "}{", u2, "}</strong>",
        "<br>Spectral density for Quantity (", ui, ") with ", xc, " coordinates:",
        "<br>  • <strong>Axis 1:</strong> ", axt1, " (", u1, ")",
        "<br>  • <strong>Axis 2:</strong> ", axt2, " (", u2, ")",
        "<br>and data (", us, "):</p>"
    )
    pretty_table(
        io,
        ustrip.(x.data),
        backend=:html,
        column_labels = x.axis2,
        row_labels=x.axis1,
        maximum_number_of_rows=10,
        vertical_crop_mode=:middle,
        style = HtmlTableStyle(;
            first_line_column_label = ["color" => "LightSeaGreen", "font-weight" => "bold"],
            row_label = ["color" => "Coral", "font-weight" => "bold"],
        )
    )
end
