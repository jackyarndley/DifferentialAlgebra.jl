struct CompiledMap{T<:Real}
    coefficients::Matrix{T}
    levels::Vector{Int}
    indices::Vector{Int}
    depth::Int
    nvars::Int
    context::Context
end
const compiledDA = CompiledMap
CompiledMap(a::Union{DA,AbstractVector{<:DA}}) = compile(a)
CompiledMap(map::CompiledMap) = copy(map)
Base.copy(map::CompiledMap) = CompiledMap(copy(map.coefficients),copy(map.levels),
    copy(map.indices),map.depth,map.nvars,map.context)
getDim(map::CompiledMap) = size(map.coefficients,1)
getOrd(map::CompiledMap) = map.depth
getVars(map::CompiledMap) = map.nvars
getTerms(map::CompiledMap) = length(map.levels)

function compile(polynomials::AbstractVector{<:DA})
    isempty(polynomials) && throw(ArgumentError("Cannot compile an empty map"))
    owners = collect(polynomials)
    ctx = valid(first(owners)); basis = ctx.basis
    all(a -> valid(a) === ctx,owners) || throw(ArgumentError("Different polynomial contexts"))
    T = mapreduce(coefftype,promote_type,owners)
    compile_map(T,owners,ctx,basis)
end
function compile_map(::Type{T},owners,ctx,basis) where T
    used = falses(last(basis.ends)); used[1] = true
    for a in owners, i in 2:a.len
        iszero(a.coeffs[i]) && continue
        k = i
        while !used[k]
            used[k] = true
            k = basis.parents[k]
        end
    end
    nodes = filter(i -> used[i],basis.traversal)
    coefficients = zeros(T,length(owners),length(nodes))
    @inbounds for (column,k) in enumerate(nodes), (row,a) in enumerate(owners)
        k <= a.len && (coefficients[row,column] = a.coeffs[k])
    end
    levels,indices = basis.degrees[nodes],basis.factors[nodes]
    CompiledMap(coefficients,levels,indices,maximum(levels),maximum(indices),ctx)
end
compile(a::DA) = compile([a])

"Evaluate a numeric map without allocations, using disjoint output and depth+1 workspace arrays."
function evaluate!(out::AbstractVector{T},map::CompiledMap,args::AbstractVector{<:Real},work::AbstractVector{T}) where {T<:Real}
    Base.require_one_based_indexing(out,args,work)
    any(x -> x isa DA,args) && throw(ArgumentError("Use evaluate for polynomial composition"))
    (Base.mightalias(out,args) || Base.mightalias(work,args) || Base.mightalias(out,work)) &&
        throw(ArgumentError("Output, inputs and workspace must not alias"))
    n = getDim(map)
    length(out) == n && length(work) >= map.depth+1 || throw(DimensionMismatch("Invalid workspace"))
    if length(args) >= map.nvars
        evaluate_numeric!(out,map,args,work,Val(true))
    else
        evaluate_numeric!(out,map,args,work,Val(false))
    end
end
function evaluate_numeric!(out::AbstractVector{T},map,args,work,::Val{COMPLETE}) where {T,COMPLETE}
    n = getDim(map)
    work[1] = one(T)
    @inbounds for j in 1:n
        out[j] = map.coefficients[j,1]
    end
    @inbounds for i in 2:length(map.levels)
        depth,var = map.levels[i],map.indices[i]
        work[depth+1] = COMPLETE || var <= length(args) ? work[depth]*args[var] : zero(T)
        for j in 1:n
            out[j] += work[depth+1]*map.coefficients[j,i]
        end
    end
    out
end
function evaluation_type(::Type{T},args) where T
    S = eltype(args)
    isconcretetype(S) ? promote_type(T,S) : foldl((R,x) -> promote_type(R,typeof(x)),args; init=T)
end
function evaluate(map::CompiledMap{T},args::AbstractVector{<:Real}) where T
    if any(x -> x isa DA,args)
        R = foldl((S,x) -> promote_type(S,x isa DA ? coefftype(x) : typeof(x)),args; init=T)
        return evaluate(map,DA{R}.(args))
    end
    R = evaluation_type(T,args)
    evaluate!(zeros(R,getDim(map)),map,args,zeros(R,map.depth+1))
end

function compose(map::CompiledMap{T},args::AbstractVector{<:DA},cutoff::Int) where T
    Base.require_one_based_indexing(args)
    ctx = map.context
    ctx.active || throw(ArgumentError("Compiled map belongs to a previous initialization"))
    all(a -> valid(a) === ctx,args) || throw(ArgumentError("Different polynomial contexts"))
    R = foldl((S,a) -> promote_type(S,coefftype(a)),args; init=T)
    if isempty(ctx.basis.products)
        compose_map(R,map,args,ctx,cutoff,Val(false))
    else
        compose_map(R,map,args,ctx,cutoff,Val(true))
    end
end
function compose_map(::Type{T},map,args,ctx,cutoff,table::Val{TABLE}) where {T,TABLE}
    n = getDim(map); capacity = ctx.basis.ends[cutoff+1]
    result = [allocate(ctx,T,capacity) for _ in 1:n]
    @inbounds for j in 1:n
        result[j].coeffs[1] = map.coefficients[j,1]
    end
    work = [allocate(ctx,T,capacity) for _ in 0:map.depth]
    work[1].coeffs[1] = one(T)
    active = [Int[] for _ in 0:map.depth]; push!(active[1],1)
    input_active = [findall(!iszero,@view a.coeffs[1:a.len]) for a in args]
    marks = zeros(Int,capacity)
    missing = allocate(ctx,T,1)
    missing_active = Int[]
    @inbounds for i in 2:length(map.levels)
        depth,var = map.levels[i],map.indices[i]
        arg = var <= length(args) ? args[var] : missing
        indices = var <= length(args) ? input_active[var] : missing_active
        sparse_multiply!(work[depth+1],active[depth+1],work[depth],active[depth],arg,indices,marks,i,cutoff,table)
        for j in 1:n
            c = map.coefficients[j,i]
            iszero(c) && continue
            for k in active[depth+1]
                value = result[j].coeffs[k]+c*work[depth+1].coeffs[k]
                result[j].coeffs[k] = keep(value,ctx) ? value : zero(T)
            end
        end
    end
    foreach(finish!,result)
    result
end

# Composition often multiplies single monomials by sparse substitutions. Track
# only touched coefficients here instead of scanning an entire dense basis.
function sparse_multiply!(out::DA{T},oi,a,ai,b,bi,marks,stamp,cutoff,::Val{TABLE}) where {T,TABLE}
    basis = out.context.basis
    @inbounds for k in oi
        out.coeffs[k] = zero(T)
    end
    empty!(oi)
    @inbounds for i in ai
        offset = basis.offsets[i]
        for j in bi
            basis.degrees[i]+basis.degrees[j] <= cutoff || continue
            k = TABLE ? Int(basis.products[offset+j]) : product_rank(basis,i,j)
            if marks[k] != stamp
                marks[k] = stamp
                push!(oi,k)
            end
            out.coeffs[k] += a.coeffs[i]*b.coeffs[j]
        end
    end
    n,last = 0,1
    @inbounds for k in oi
        if keep(out.coeffs[k],out.context)
            n += 1; oi[n] = k; last = max(last,k)
        else
            out.coeffs[k] = zero(T)
        end
    end
    resize!(oi,n)
    out.len = last
    out
end
evaluate(map::CompiledMap,args::AbstractVector{<:DA}) = compose(map,args,map.context.cutoff)
evaluate(a::DA,args::AbstractVector{<:Real}) = only(evaluate(compile(a),args))
evaluate(a::AbstractVector{<:DA},args::AbstractVector{<:Real}) = evaluate(compile(a),args)
evalScalar(a::DA,value::Real) = evaluate(a,[value])
evalScalar(a::Union{CompiledMap,AbstractVector{<:DA}},value::Real) = evaluate(a,[value])
function evaluate(map::CompiledMap,args::AbstractVector{<:Real},result::AbstractVector)
    values = evaluate(map,args)
    length(values) == length(result) || throw(DimensionMismatch("Invalid output length"))
    copyto!(result,values)
end

function _linear_transform(A::AbstractMatrix{<:Real},x::AbstractVector{<:DA})
    size(A,2) == length(x) || throw(DimensionMismatch("Linear map dimensions differ"))
    ctx = valid(first(x))
    T = foldl((R,a) -> promote_type(R,coefftype(a)),x; init=eltype(A))
    result = [allocate(ctx,T) for _ in axes(A,1)]
    for j in eachindex(x), i in eachindex(result)
        c = A[i,j]
        iszero(c) && continue
        weighted_sum!(result[i],result[i],one(T),x[j],c)
    end
    result
end
function invert(f::AbstractVector{<:DA})
    isempty(f) && throw(DimensionMismatch("An empty map cannot be inverted"))
    ctx = valid(first(f)); nv = ctx.basis.variables
    1 <= length(f) <= nv || throw(DimensionMismatch("Map dimension exceeds the independent variables"))
    all(a -> valid(a) === ctx,f) || throw(ArgumentError("Different polynomial contexts"))
    T = mapreduce(coefftype,promote_type,f)
    if length(f) < nv
        return invert(vcat(f,[variable(i,T) for i in length(f)+1:nv]))[1:length(f)]
    end
    x = [variable(i,T) for i in 1:nv]
    c,A = cons.(f),linear(f)
    Ainv = lu(A) \ Matrix{T}(I,nv,nv)
    linear_inverse = _linear_transform(Ainv,x)
    nonlinear = compile(_linear_transform(Ainv,trim.(f,2)))
    result = linear_inverse
    # Explicit local orders avoid changing the context or racing other readers.
    for degree in 2:ctx.cutoff
        result = linear_inverse - compose(nonlinear,result,degree)
    end
    compose(compile(result),x-c,ctx.cutoff)
end
