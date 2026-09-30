"DACE text representation; BigFloat coefficients retain all their decimal digits."
function toString(a::DA)
    valid(a)
    io = IOBuffer()
    if iszero(a)
        println(io,"        ALL COEFFICIENTS ZERO")
    else
        println(io,"     I  COEFFICIENT              ORDER EXPONENTS")
        for (i,m) in enumerate(getMonomials(a))
            if m.coefficient isa Union{Float32,Float64}
                @printf(io,"%6d  %24.16e%4d ",i,m.coefficient,order(m))
            else
                print(io,i,"  ",m.coefficient,"  ",order(m)," ")
            end
            for power in m.exponents; @printf(io," %2d",power); end
            println(io)
        end
    end
    println(io,"------------------------------------------------")
    String(take!(io))
end

"""
    fromString(text, T=Float64)

Read a DACE text polynomial into the current initialization. Terms above the
initialized order, or involving extra variables, are discarded as in the DACE C core.
Use `T=BigFloat` inside `setprecision` to retain arbitrary-precision coefficients.
"""
fromString(text::AbstractString,::Type{T}=Float64) where {T<:Real} = fromString(split(text,'\n'),T)
function fromString(lines::AbstractVector{<:AbstractString},::Type{T}=Float64) where {T<:Real}
    nonempty = filter(!isempty,strip.(lines))
    isempty(nonempty) && throw(ArgumentError("Empty polynomial text"))
    a = DA{T}()
    occursin(r"^ALL (COEFFICIENTS|COMPONENTS) ZERO$",nonempty[1]) && return a
    occursin(r"^I\s+COEFFICIENT\s+ORDER EXPONENTS$",nonempty[1]) || throw(ArgumentError("Invalid DACE text header"))
    for line in nonempty[2:end]
        startswith(line,"---") && break
        fields = split(line)
        length(fields) >= 4 || throw(ArgumentError("Invalid DACE coefficient line"))
        index = parse(Int,fields[1])
        value = parse(T,replace(fields[2],'D'=>'E','d'=>'e'))
        degree = parse(Int,fields[3])
        powers = parse.(Int,fields[4:end])
        index > 0 && all(>=(0),powers) && sum(big,powers) == degree || throw(ArgumentError("Invalid DACE monomial"))
        b = a.context.basis
        degree <= b.order && all(iszero,@view powers[min(end,b.variables)+1:end]) || continue
        setCoefficient!(a,powers,value)
    end
    a
end
Base.parse(::Type{DA},text::AbstractString) = fromString(text)
Base.parse(::Type{DA{T}},text::AbstractString) where {T<:Real} = fromString(text,T)

# Interoperate with the packed C-core blob: five UInt32 header words followed
# by (UInt32,UInt32,Float64) monomials. Zero polynomials include one dummy term.
# Existing DACE blobs use native endianness; accept either byte order on input
# and write the little-endian format used by x86 Linux and Windows releases.
const BLOB_MAGIC = UInt32(0x1e304144)
function blob_index(powers,order)
    index = big(0)
    for power in Iterators.reverse(powers)
        index = index*(order+1)+power
    end
    index <= typemax(UInt32) || throw(ArgumentError("This basis exceeds the legacy DACE blob index range; use text IO"))
    UInt32(index)
end
function Base.write(io::IO,a::DA{Float64})
    b = valid(a).basis
    nv1 = cld(b.variables,2)
    terms = getMonomials(a)
    # Validate before writing anything to the caller's stream.
    indices = [(blob_index(@view(m.exponents[1:nv1]),b.order),
                blob_index(@view(m.exponents[nv1+1:end]),b.order)) for m in terms]
    bytes = 0
    for value in (BLOB_MAGIC,UInt32(b.order),UInt32(nv1),UInt32(b.variables-nv1),UInt32(length(terms)))
        bytes += write(io,htol(value))
    end
    for (m,(i,j)) in zip(terms,indices)
        bytes += write(io,htol(i),htol(j),htol(reinterpret(UInt64,m.coefficient)))
    end
    isempty(terms) && (bytes += write(io,zeros(UInt8,16)))
    bytes
end
function Base.write(io::IO,a::DA)
    throw(ArgumentError("Legacy binary DACE stores Float64; use toString/fromString for other coefficient types"))
end
function Base.read(io::IO,::Type{DA{T}}) where {T<:Real}
    marker = read(io,UInt32)
    decode = marker == ltoh(BLOB_MAGIC) ? ltoh : marker == ntoh(BLOB_MAGIC) ? ntoh :
        throw(ArgumentError("Invalid DACE blob magic"))
    order,nv1,nv2,count = (Int(decode(read(io,UInt32))) for _ in 1:4)
    1 <= order <= 65535 && 1 <= nv1+nv2 <= 1024 || throw(ArgumentError("Invalid DACE blob dimensions"))
    count <= binomial(big(order)+nv1+nv2,nv1+nv2) || throw(ArgumentError("Invalid DACE blob term count"))
    a = DA{T}()
    powers = zeros(Int,nv1+nv2)
    for _ in 1:max(count,1)
        i,j = decode(read(io,UInt32)),decode(read(io,UInt32))
        value = reinterpret(Float64,decode(read(io,UInt64)))
        count == 0 && continue
        for (range,index) in ((1:nv1,i),(nv1+1:nv1+nv2,j))
            for k in range
                index,powers[k] = divrem(index,order+1)
            end
            index == 0 || throw(ArgumentError("Invalid DACE blob exponent index"))
        end
        sum(powers) <= order || throw(ArgumentError("DACE blob monomial exceeds its source order"))
        b = a.context.basis
        sum(powers) <= b.order && all(iszero,@view powers[min(end,b.variables)+1:end]) || continue
        setCoefficient!(a,powers,value)
    end
    a
end
Base.read(io::IO,::Type{DA}) = read(io,DA{Float64})
