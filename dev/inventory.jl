using SHA

# Parse source only: an inventory must not load optional dependencies or run examples.
const SOURCE_DIRECTORIES = ("src", "ext", "examples", "test", "docs", "benchmark", "dev")
const INVENTORY_COLUMNS = ("file", "line", "kind", "name", "signature", "body_sha1", "context")

function source_files(root)
    files = String[]
    for directory in SOURCE_DIRECTORIES
        path = joinpath(root, directory)
        isdir(path) || continue
        for (folder, children, names) in walkdir(path)
            filter!(name -> name ∉ ("build", "generated", "output", ".git"), children)
            append!(files, (joinpath(folder, name) for name in names if endswith(name, ".jl")))
        end
    end
    return sort!(files)
end

function without_lines(value)
    value isa Expr || return value
    # Macro calls reserve their second argument for a source location. Removing
    # that slot would change the syntax and could hide the macro's first argument.
    if value.head === :macrocall
        return Expr(
            :macrocall, first(value.args), LineNumberNode(0),
            (without_lines(arg) for arg in value.args[3:end])...
        )
    end
    return Expr(value.head, (without_lines(arg) for arg in value.args if !(arg isa LineNumberNode))...)
end

source_text(value) = replace(string(without_lines(value)), r"\s+" => " ")

function signature_call(signature)
    while signature isa Expr && signature.head in (:where, :(::))
        signature = first(signature.args)
    end
    return signature isa Expr && signature.head === :call ? signature : nothing
end

function function_name(signature)
    call = signature_call(signature)
    call === nothing && return source_text(signature)
    name = first(call.args)
    return name isa Expr && name.head === :(::) ? "call " * source_text(name) : source_text(name)
end

function first_line(body, fallback)
    body isa Expr && body.head === :block || return fallback
    for arg in body.args
        arg isa LineNumberNode && return Int(arg.line)
    end
    return fallback
end

function inventory_entry!(rows, file, line, kind, signature, body, context, template)
    label = template ? "template-" * kind : kind
    name = kind == "anonymous" ? "(anonymous)" : function_name(signature)
    digest = bytes2hex(sha1(repr(without_lines(body))))
    push!(
        rows, (
            file, string(first_line(body, line)), label, name,
            source_text(signature), digest, isempty(context) ? "(top-level)" : join(context, " / "),
        )
    )
    return name
end

function inventory!(rows, node, file; line = 1, context = String[], template = false)
    node isa Expr || return nothing
    if node.head === :function
        signature = first(node.args)
        body = length(node.args) == 1 ? nothing : last(node.args)
        kind = body === nothing ? "declaration" :
            signature isa Expr && signature.head === :tuple ? "anonymous" : "method"
        name = inventory_entry!(rows, file, line, kind, signature, body, context, template)
        inventory!(rows, signature, file; line, context = [context; name], template)
        inventory!(rows, body, file; line, context = [context; name], template)
        return nothing
    elseif node.head === :(=) && signature_call(first(node.args)) !== nothing
        signature, body = node.args
        name = inventory_entry!(rows, file, line, "method", signature, body, context, template)
        inventory!(rows, signature, file; line, context = [context; name], template)
        inventory!(rows, body, file; line, context = [context; name], template)
        return nothing
    elseif node.head === :->
        signature, body = node.args
        inventory_entry!(rows, file, line, "anonymous", signature, body, context, template)
        inventory!(rows, body, file; line, context = [context; "(anonymous)"], template)
        return nothing
    elseif node.head === :macro
        signature, body = node.args
        name = inventory_entry!(rows, file, line, "macro", signature, body, context, template)
        inventory!(rows, body, file; line, context = [context; name], template)
        return nothing
    elseif node.head === :macrocall
        macro_name = first(node.args)
        if macro_name isa Expr && macro_name.head === :.
            macro_name = last(macro_name.args)
        end
        macro_name isa QuoteNode && (macro_name = macro_name.value)
        template |= macro_name === Symbol("@eval")
    elseif node.head === :for
        context = [context; "for " * source_text(first(node.args))]
    elseif node.head === :module
        context = [context; "module " * source_text(node.args[2])]
    elseif node.head === :quote && !template
        # Quoted code is data until explicitly evaluated. Report its definitions
        # as templates rather than claiming that the methods exist at runtime.
        template = true
    elseif node.head in (:error, :incomplete)
        throw(ArgumentError("Cannot inventory invalid Julia source in $file at line $line"))
    end
    for arg in node.args
        if arg isa LineNumberNode
            line = Int(arg.line)
        else
            inventory!(rows, arg, file; line, context, template)
        end
    end
    return nothing
end

function function_inventory(root)
    rows = NTuple{7, String}[]
    files = source_files(root)
    for file in files
        parsed = Meta.parseall(read(file, String); filename = file)
        inventory!(rows, parsed, replace(relpath(file, root), '\\' => '/'))
    end
    sort!(rows; by = row -> (row[1], parse(Int, row[2]), row[3], row[5]))
    return files, rows
end

function write_inventory(io, rows)
    println(io, join(INVENTORY_COLUMNS, '\t'))
    for row in rows
        println(io, join(row, '\t'))
    end
    return nothing
end

function main(args)
    length(args) <= 1 || throw(ArgumentError("Usage: julia dev/inventory.jl [output.tsv]"))
    root = dirname(@__DIR__)
    files, rows = function_inventory(root)
    if isempty(args)
        write_inventory(stdout, rows)
    else
        open(io -> write_inventory(io, rows), only(args), "w")
    end
    println(stderr, "Inventoried ", length(rows), " definitions in ", length(files), " Julia files.")
    for kind in sort!(unique(row[3] for row in rows))
        println(stderr, "  ", kind, ": ", count(row -> row[3] == kind, rows))
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
