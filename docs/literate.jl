using Base64: base64decode

# Makie's rich HTML display embeds PNG data. Convert that display to a normal
# Markdown image so Documenter resolves it with both flat and pretty URLs.
function example_images(markdown, destination, name)
    index = 0
    pattern = r"```@raw html\s*<img\b[^>]*\bsrc=\"data:image/png;base64,\s*([A-Za-z0-9+/=\s]+)\"\s*/?>\s*```"
    return replace(
        markdown, pattern => function (block)
            encoded = only(match(pattern, block).captures)
            index += 1
            filename = "$name-figure-$index.png"
            write(joinpath(destination, filename), base64decode(encoded))
            return "![]($filename)"
        end
    )
end
