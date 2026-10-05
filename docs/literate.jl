using Base64: base64decode

# The overview supplies both the execution order and topic groups in navigation.
function example_groups(markdown)
    groups = Pair{String, Vector{Pair{String, String}}}[]
    topic = ""
    for line in split(markdown, '\n')
        heading = match(r"^## (.+)$", strip(line))
        if heading !== nothing
            topic = String(heading.captures[1])
            continue
        end
        entry = match(r"\[([^\]]+)\]\(generated/([^.)]+)\.md\)", line)
        entry === nothing && continue
        isempty(topic) && error("Each tutorial must belong to a topic heading")
        if isempty(groups) || first(last(groups)) != topic
            push!(groups, topic => Pair{String, String}[])
        end
        push!(last(last(groups)), String(entry.captures[1]) => String(entry.captures[2]))
    end
    return groups
end

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
