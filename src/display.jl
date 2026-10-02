# Human-readable expressions. Display is deliberately independent of file formats.
# Use scalar show methods so coefficient precision and custom types are preserved.
function show_term(io::IO, coefficient::Real, exponents, names, first::Bool)
    negative = coefficient < zero(coefficient)
    if first
        negative && print(io, '-')
    else
        print(io, negative ? " - " : " + ")
    end
    # Widen negative integers before negation, including typemin(Int).
    magnitude = negative ? -(coefficient isa Integer ? big(coefficient) : coefficient) : coefficient
    show(io, magnitude)
    any(!iszero, exponents) && print(io, ' ')
    for (i, power) in enumerate(exponents)
        iszero(power) && continue
        print(io, names[i])
        if power != 1
            for digit in string(power)
                print(io, ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')[digit - '0' + 1])
            end
        end
    end
    return nothing
end

function Base.show(io::IO, p::TaylorPolynomial)
    if !p.algebra.active
        print(io, "TaylorPolynomial(inactive)")
        return nothing
    end
    basis = p.algebra.basis
    # Honor Julia's :limit context in the REPL and array displays. Ordinary
    # show/string calls retain every term; compact mode still shows the variables.
    limit = get(io, :limit, false) ? (get(io, :compact, false) ? 4 : 8) : typemax(Int)
    terms = 0
    for i in 1:p.len
        iszero(p.coeffs[i]) && continue
        if terms == limit
            print(io, " + …")
            return nothing
        end
        show_term(io, p.coeffs[i], @view(basis.exponents[:, i]), p.algebra.names, terms == 0)
        terms += 1
    end
    terms == 0 && show(io, zero(coefficient_type(p)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", p::TaylorPolynomial)
    p.algebra.active || return show(io, p)
    n = p.algebra.basis.variables
    print(
        io, typeof(p), " polynomial in ", n, n == 1 ? " variable" : " variables",
        " (order ≤ ", p.algebra.basis.order, "):\n  "
    )
    return show(io, p)
end

function Base.show(io::IO, m::Monomial)
    print(io, "Monomial(")
    show(io, m.coefficient)
    print(io, ", [")
    join(io, m.exponents, ", ")
    return print(io, "])")
end
