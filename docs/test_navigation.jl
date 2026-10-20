using Test

# Check destinations and referenced assets without depending on Documenter's
# sidebar markup or requiring every standalone tutorial to produce a figure.
function test_navigation(examples, buildroot, generated)
    buildroot = abspath(buildroot)
    pagepath(name) = isfile(joinpath(buildroot, name * ".html")) ?
        joinpath(buildroot, name * ".html") : joinpath(buildroot, name, "index.html")
    @testset "Built catalog, local links and assets" begin
        for (_, name) in examples
            page = pagepath("generated/" * name)
            @test isfile(page)
            html = read(page, String)
            sources = [only(m.captures) for m in eachmatch(r"<img\b[^>]*src=\"([^\"]+)\"", html)]
            markdown = read(joinpath(generated, name * ".md"), String)
            for image in eachmatch(r"!\[[^\]]*\]\(([^)]+\.png)\)", markdown)
                filename = only(image.captures)
                @test isfile(joinpath(generated, filename))
                @test any(href -> endswith(href, filename) && isfile(joinpath(dirname(page), href)), sources)
            end
        end
        pages = [
            joinpath(dir, file) for (dir, _, files) in walkdir(buildroot)
                for file in files if endswith(file, ".html")
        ]
        @test !isempty(pages)
        for page in pages
            html = read(page, String)
            for link in eachmatch(r"\b(?:href|src)=\"([^\"]*)\"", html)
                href = only(link.captures)
                # External URLs, embedded displays and fragment-only links do
                # not name a local file. Documenter checks internal references.
                (
                    isempty(href) || startswith(href, '#') || startswith(href, "//") ||
                        occursin(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", href)
                ) && continue
                path = first(split(href, r"[?#]"))
                target = startswith(path, '/') ? joinpath(buildroot, path[2:end]) :
                    normpath(joinpath(dirname(page), path))
                @test first(splitpath(relpath(target, buildroot))) != ".."
                @test isfile(target) || isfile(joinpath(target, "index.html"))
            end
        end
    end
    return nothing
end
