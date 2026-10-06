using Test
isdefined(@__MODULE__, :example_groups) || include("literate.jl")

# Check rendered pages, including direct links to a tutorial in another topic.
# Catalog registration alone cannot detect a collapsed or inconsistent sidebar.
@testset "Unified built documentation navigation and figures" begin
    groups = example_groups(read(joinpath(@__DIR__, "src/examples.md"), String))
    entries = reduce(vcat, last.(groups))
    buildroot = joinpath(@__DIR__, "build")
    pagepath(name) = isfile(joinpath(buildroot, name * ".html")) ? joinpath(buildroot, name * ".html") : joinpath(buildroot, name, "index.html")
    resolve(page, href) = begin
        target = isempty(href) ? page : normpath(joinpath(dirname(page), first(split(href, '#'))))
        endswith(target, ".html") ? target : joinpath(target, "index.html")
    end
    pages = [joinpath(dir, file) for (dir, _, files) in walkdir(buildroot) for file in files if endswith(file, ".html") && !startswith(relpath(dir, buildroot), "assets")]
    @test length(pages) > length(entries)
    for page in pages
        html = read(page, String)
        sidebar = match(r"<nav class=\"docs-sidebar\">(.*?)</nav>"s, html)
        @test sidebar !== nothing
        sidebar === nothing && continue
        menu = only(sidebar.captures)
        # At collapselevel=3, topic groups have always-visible lists. Only
        # an individual page's longer contents may use a collapsed toggle.
        expanded = all(topic -> occursin("<span class=\"tocitem\">$topic</span><ul>", menu), first.(groups))
        @test expanded
        # Documenter uses a valueless href for the active page with pretty URLs.
        links = Dict(m.captures[2] => something(m.captures[1], "") for m in eachmatch(r"<a class=\"tocitem\" href(?:=\"([^\"]*)\")?>([^<]+)</a>", menu))
        registered = all(entry -> haskey(links, first(entry)), entries)
        resolved = all(entry -> haskey(links, first(entry)) && resolve(page, links[first(entry)]) == normpath(pagepath("generated/" * last(entry))), entries)
        @test registered
        @test resolved
    end
    for (_, name) in entries
        markdown = read(joinpath(@__DIR__, "src/generated", name * ".md"), String)
        images = collect(eachmatch(r"!\[\]\(([^)]+\.png)\)", markdown))
        @test !isempty(images)
        page = pagepath("generated/" * name)
        html = read(page, String)
        links = [only(m.captures) for m in eachmatch(r"<img\b[^>]*src=\"([^\"]+)\"", html)]
        @test all(images) do image
            filename = only(image.captures)
            any(links) do href
                endswith(href, filename) && isfile(normpath(joinpath(dirname(page), href)))
            end
        end
    end
end
