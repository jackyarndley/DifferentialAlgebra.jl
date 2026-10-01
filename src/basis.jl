"Degree-ordered monomials and a bounded, packed multiplication lookup table."
struct MonomialBasis
    order::Int
    variables::Int
    table_bytes::Int
    exponents::Matrix{Int}
    degrees::Vector{Int}
    ends::Vector{Int}
    counts::Matrix{Int}
    offsets::Vector{Int}
    products::Vector{Int32}
    splits::Vector{Int}
    left::Vector{Int32}
    right::Vector{Int32}
    parents::Vector{Int}
    factors::Vector{Int}
    traversal::Vector{Int}
end

# Recurrence pairs are sorted by the left monomial. Trim the range once so
# coefficient reductions have no early exit and can use SIMD for machine floats.
@inline function split_end(b::MonomialBasis, k::Int, len::Int)
    hi = b.splits[k + 1]
    len >= k && return hi
    lo = b.splits[k] + 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if b.left[mid] <= len
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return hi
end

# Rank a reverse-lexicographic monomial without allocating a tuple or dictionary.
@inline function rank(b::MonomialBasis, exponents, degree::Int)
    index = degree == 0 ? 1 : b.ends[degree] + 1
    remaining = degree
    @inbounds for v in 1:(b.variables - 1)
        remaining -= Int(exponents[v])
        remaining > 0 && (index += b.counts[remaining, b.variables - v + 1])
    end
    return index
end
@inline function product_rank(b::MonomialBasis, i::Int, j::Int)
    degree = b.degrees[i] + b.degrees[j]
    index = degree == 0 ? 1 : b.ends[degree] + 1
    @inbounds for v in 1:(b.variables - 1)
        degree -= b.exponents[v, i] + b.exponents[v, j]
        degree > 0 && (index += b.counts[degree, b.variables - v + 1])
    end
    return index
end

function foreach_monomial(f, no::Int, nv::Int)
    current = zeros(Int, nv)
    function visit(v, remaining, degree)
        return if v == nv
            current[v] = remaining
            f(current, degree)
        else
            for value in remaining:-1:0
                current[v] = value
                visit(v + 1, remaining - value, degree)
            end
        end
    end
    for degree in 0:no
        visit(1, degree, degree)
    end
    return nothing
end

function MonomialBasis(no::Int, nv::Int; table_bytes::Int = 32 * 1024^2)
    ends = [Int(binomial(big(nv + d), nv)) for d in 0:no]
    nm = last(ends)
    counts = [Int(binomial(big(d + v), v)) for d in 0:no, v in 0:nv]
    exponents = Matrix{Int}(undef, nv, nm)
    degrees = Vector{Int}(undef, nm)
    column = 0
    foreach_monomial(no, nv) do powers, degree
        column += 1
        exponents[:, column] = powers
        degrees[column] = degree
    end
    current = zeros(Int, nv)
    offsets, products = zeros(Int, nm), Int32[]
    parents, factors = ones(Int, nm), zeros(Int, nm)
    traversal = Int[]
    splits, left, right = zeros(Int, nm + 1), Int32[], Int32[]
    b = MonomialBasis(no, nv, table_bytes, exponents, degrees, ends, counts, offsets, products, splits, left, right, parents, factors, traversal)
    pairs = sum(i -> ends[no - degrees[i] + 1], 1:nm)
    if pairs <= table_bytes ÷ (3sizeof(Int32))
        sizehint!(products, pairs)
        for i in 1:nm
            offsets[i] = length(products)
            for j in 1:ends[no - degrees[i] + 1]
                push!(products, Int32(product_rank(b, i, j)))
            end
        end
        # Group nonconstant-left products by output monomial for recurrences.
        for i in 2:nm, j in 1:ends[no - degrees[i] + 1]
            k = Int(products[offsets[i] + j])
            splits[k + 1] += 1
        end
        cumsum!(splits, splits)
        resize!(left, last(splits)); resize!(right, last(splits))
        positions = copy(splits)
        for i in 2:nm, j in 1:ends[no - degrees[i] + 1]
            k = Int(products[offsets[i] + j])
            positions[k] += 1
            left[positions[k]], right[positions[k]] = i, j
        end
    end
    children = [Int[] for _ in 1:nm]
    for i in 2:nm
        current .= @view exponents[:, i]
        v = findlast(!iszero, current)::Int
        current[v] -= 1
        parents[i], factors[i] = rank(b, current, degrees[i] - 1), v
        push!(children[parents[i]], i)
    end
    stack = [1]
    while !isempty(stack)
        i = pop!(stack)
        push!(traversal, i)
        append!(stack, Iterators.reverse(children[i]))
    end
    return b
end
