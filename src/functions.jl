# Compose univariate Taylor coefficients with the nonconstant part of a DA.
# Horner stages need only successively increasing total degrees, and reuse buffers.
function series(a::DA,coefficients::AbstractVector{T}) where T
    ctx = valid(a)
    R = promote_type(coefftype(a),T)
    a.len == 1 && return DA{R}(coefficients[1])
    n = min(ctx.cutoff,length(coefficients)-1)
    h = DA{R}(a); h.coeffs[1] = zero(R)
    result, scratch = allocate(ctx,R), allocate(ctx,R)
    result.coeffs[1] = coefficients[n+1]
    for k in n-1:-1:0
        multiply!(scratch,h,result,n-k)
        scratch.coeffs[1] += coefficients[k+1]
        finish!(scratch,scratch.len)
        result,scratch = scratch,result
    end
    result
end
function domain_call(f,x)
    try
        f(x)
    catch err
        err isa DomainError || rethrow()
        throw(DAError(sprint(showerror,err)))
    end
end

# Euler's homogeneous derivative E(x^alpha)=|alpha|x^alpha turns
# E(exp(a))=exp(a)E(a), a E(log(a))=E(a), and a E(a^p)=p a^p E(a)
# into triangular coefficient recurrences. One traversal replaces repeated
# polynomial products, with no cutoff applied to intermediate Taylor factors.
function recurrence(a::DA,c0,p,::Val{MODE}) where MODE
    ctx = valid(a); basis = ctx.basis
    T = promote_type(coefftype(a),typeof(c0),typeof(p))
    a.len == 1 && return DA{T}(c0)
    n = basis.ends[ctx.cutoff+1]
    out = allocate_undef(ctx,T,n); out.coeffs[1] = c0
    a0 = a.coeffs[1]
    @inbounds for k in 2:n
        degree = basis.degrees[k]
        value = zero(T)
        for t in basis.splits[k]+1:basis.splits[k+1]
            i,j = Int(basis.left[t]),Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            iszero(ai) && continue
            weight = MODE === :exp ? basis.degrees[i] : MODE === :log ? -basis.degrees[j] : p*basis.degrees[i]-basis.degrees[j]
            value += weight*ai*out.coeffs[j]
        end
        if MODE === :log
            value += k <= a.len ? degree*a.coeffs[k] : zero(T)
        end
        out.coeffs[k] = MODE === :exp ? value/degree : value/(degree*a0)
    end
    finish!(out,n)
end

function trig_recurrence(a::DA,hyperbolic::Bool)
    x = cons(a); ctx = a.context; basis = ctx.basis
    s0,c0 = hyperbolic ? (sinh(x),cosh(x)) : sincos(x)
    T = typeof(s0)
    a.len == 1 && return DA{T}(s0),DA{T}(c0)
    n = basis.ends[ctx.cutoff+1]
    s,c = allocate_undef(ctx,T,n),allocate_undef(ctx,T,n)
    s.coeffs[1],c.coeffs[1] = s0,c0
    @inbounds for k in 2:n
        sv,cv = zero(T),zero(T)
        for t in basis.splits[k]+1:basis.splits[k+1]
            i,j = Int(basis.left[t]),Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            iszero(ai) && continue
            weight = basis.degrees[i]*ai
            sv += weight*c.coeffs[j]
            cv += weight*s.coeffs[j]
        end
        s.coeffs[k] = sv/basis.degrees[k]
        c.coeffs[k] = (hyperbolic ? cv : -cv)/basis.degrees[k]
    end
    finish!(s,n),finish!(c,n)
end

function Base.exp(a::DA)
    x = cons(a); c0 = exp(x)
    !isempty(a.context.basis.products) && return recurrence(a,c0,zero(c0),Val(:exp))
    c = Vector{typeof(c0)}(undef,a.context.cutoff+1); c[1] = c0
    @inbounds for k in 1:length(c)-1
        c[k+1] = c[k]/k
    end
    series(a,c)
end
function Base.log(a::DA)
    x = cons(a); c0 = domain_call(log,x)
    iszero(x) && throw(DAError("Logarithm requires a nonzero constant"))
    !isempty(a.context.basis.products) && return recurrence(a,c0,zero(c0),Val(:log))
    c = Vector{typeof(c0)}(undef,a.context.cutoff+1); c[1] = c0
    factor = one(x)
    @inbounds for k in 1:length(c)-1
        c[k+1] = factor/k
        factor = -factor
    end
    series(a/x,c)
end
Base.log2(a::DA) = log(a)/log(convert(typeof(float(cons(a))),2))
Base.log10(a::DA) = log(a)/log(convert(typeof(float(cons(a))),10))
Base.log(b::Real,a::DA) = log(a)/log(b)
Base.log(::Irrational{:ℯ},a::DA) = log(a)

function power_series(a::DA,p::Real,c0)
    x = cons(a)
    iszero(x) && throw(DAError("Noninteger powers require a nonzero expansion point"))
    !isempty(a.context.basis.products) && return recurrence(a,c0,p,Val(:power))
    T = promote_type(typeof(c0),typeof(p))
    c = Vector{T}(undef,a.context.cutoff+1); c[1] = c0
    @inbounds for k in 1:length(c)-1
        c[k+1] = c[k]*(p-(k-1))/k
    end
    series(a/x,c)
end
function Base.:^(a::DA,p::AbstractFloat)
    R = promote_type(coefftype(a),typeof(p))
    isinteger(p) && typemin(Int32) < p <= typemax(Int32) && return convert(DA{R},a)^Int(p)
    iszero(a) && p > 0 && return zero(DA{R})
    power_series(a,p,domain_call(x -> x^p,cons(a)))
end
Base.:^(a::DA,p::Rational) = a^convert(float(promote_type(coefftype(a),typeof(p))),p)
powd(a::DA,p::Real) = a^float(p)
Base.:^(a::Real,b::DA) = exp(log(a)*b)
Base.:^(a::DA,b::DA) = exp(log(a)*b)
Base.:^(::Irrational{:ℯ},b::DA) = exp(b)
@inline coefficient_at(a::DA,k::Int) = k <= a.len ? a.coeffs[k] : zero(coefftype(a))
@inline coefficient_at(a::Real,k::Int) = k == 1 ? a : zero(a)

# b*y=a is triangular by total degree. Division needs one recurrence, avoiding
# both the homogeneous-degree weights of a generic power and a second product.
function quotient(a::Real,b::DA)
    ctx = a isa DA ? compatible(a,b) : valid(b)
    b0 = cons(b)
    iszero(b0) && throw(DAError("Division requires a nonzero denominator constant"))
    if b.len == 1
        result = a/b0
        return result isa DA ? result : DA{typeof(result)}(result)
    end
    basis = ctx.basis
    if isempty(basis.products)
        return a*power_series(b,-one(b0),inv(b0))
    end
    T = promote_type(a isa DA ? coefftype(a) : typeof(a),coefftype(b))
    T = typeof(one(T)/one(T))
    n = basis.ends[ctx.cutoff+1]
    result = allocate_undef(ctx,T,n)
    result.coeffs[1] = cons(a)/b0
    @inbounds for k in 2:n
        value = convert(T,coefficient_at(a,k))
        for t in basis.splits[k]+1:basis.splits[k+1]
            i,j = Int(basis.left[t]),Int(basis.right[t])
            i > b.len && break
            bi = b.coeffs[i]
            iszero(bi) && continue
            value = muladd(-bi,result.coeffs[j],value)
        end
        result.coeffs[k] = value/b0
    end
    finish!(result,n)
end
Base.inv(a::DA) = quotient(one(coefftype(a)),a)

# y^2=a: use each unordered pair of nonconstant monomials only once.
function square_root(a::DA,c0)
    ctx = valid(a); basis = ctx.basis
    a.len == 1 && return DA(c0)
    n = basis.ends[ctx.cutoff+1]
    result = allocate_undef(ctx,typeof(c0),n)
    result.coeffs[1] = c0
    denominator = c0+c0
    @inbounds for k in 2:n
        value = coefficient_at(a,k)
        for t in basis.splits[k]+1:basis.splits[k+1]
            i,j = Int(basis.left[t]),Int(basis.right[t])
            i > j && break
            term = result.coeffs[i]*result.coeffs[j]
            value -= i == j ? term : term+term
        end
        result.coeffs[k] = value/denominator
    end
    finish!(result,n)
end
function root(a::DA,p::Integer=2)
    p != 0 || throw(DomainError(p,"Zeroth root is undefined"))
    p == 1 && return copy(a)
    p > 0 && iszero(a) && return zero(a)
    x = cons(a)
    iszero(x) && throw(DAError("Root is not analytic at a zero constant"))
    iseven(p) && x < 0 && throw(DAError("Even root of a negative constant"))
    exponent = one(float(x))/p
    c0 = p == 2 ? sqrt(x) : p == 3 ? cbrt(x) : copysign(abs(x)^exponent,x)
    p == 2 && !isempty(a.context.basis.products) && return square_root(a,c0)
    power_series(a,exponent,c0)
end
Base.sqrt(a::DA) = root(a,2)
Base.cbrt(a::DA) = root(a,3)
isrt(a::DA) = root(a,-2)
icrt(a::DA) = root(a,-3)
function root(x::Real,p::Integer=2)
    p != 0 || throw(DomainError(p,"Zeroth root is undefined"))
    p == 1 && return x
    exponent = one(float(x))/p
    isodd(p) ? copysign(abs(x)^exponent,x) : x^exponent
end
isrt(x::Real) = inv(sqrt(x))
icrt(x::Real) = inv(cbrt(x))

for fn in (:sin,:cos,:sinh,:cosh)
    @eval function Base.$fn(a::DA)
        if !isempty(valid(a).basis.products)
            pair = trig_recurrence(a,$(fn in (:sinh,:cosh)))
            return pair[$(fn in (:sin,:sinh) ? 1 : 2)]
        end
        x = cons(a)
        s,c = $(fn in (:sin,:cos) ? :(sincos(x)) : :((sinh(x),cosh(x))))
        values = $(fn == :sin ? :((s,c,-s,-c)) : fn == :cos ? :((c,-s,-c,s)) : fn == :sinh ? :((s,c,s,c)) : :((c,s,c,s)))
        coeffs = Vector{typeof(s)}(undef,a.context.cutoff+1)
        factor = one(s)
        @inbounds for k in 0:length(coeffs)-1
            coeffs[k+1] = values[mod(k,4)+1]*factor
            factor /= k+1
        end
        series(a,coeffs)
    end
end
Base.sincos(a::DA) = isempty(valid(a).basis.products) ? (sin(a),cos(a)) : trig_recurrence(a,false)
for (fn,sign) in ((:tan,1),(:tanh,-1))
    @eval function Base.$fn(a::DA)
        x = cons(a); c0 = $fn(x)
        c = zeros(typeof(c0),a.context.cutoff+1); c[1] = c0
        @inbounds for k in 0:length(c)-2
            value = k == 0 ? one(c0) : zero(c0)
            for j in 0:k
                value += $sign*c[j+1]*c[k-j+1]
            end
            c[k+2] = value/(k+1)
        end
        series(a,c)
    end
end

# If y=(q0+q1*t+q2*t^2)^p, q*y'=p*q'*y gives this linear recurrence.
function quadratic_power(q0,q1,q2,p,n)
    iszero(q0) && throw(DAError("Singular derivative at the expansion point"))
    c0 = domain_call(x -> x^p,q0)
    c = zeros(typeof(c0),n+1); c[1] = c0
    @inbounds for k in 0:n-1
        value = (p-k)*q1*c[k+1]
        k > 0 && (value += (2p-k+1)*q2*c[k])
        c[k+2] = value/((k+1)*q0)
    end
    c
end
function integral_series(a::DA,c0,derivative)
    T = promote_type(typeof(c0),eltype(derivative))
    c = Vector{T}(undef,length(derivative)+1); c[1] = c0
    @inbounds for k in eachindex(derivative)
        c[k+1] = derivative[k]/k
    end
    series(a,c)
end
for fn in (:asin,:acos,:atan,:asinh,:acosh,:atanh)
    @eval function Base.$fn(a::DA)
        x = cons(a); u = one(float(x)); c0 = domain_call($fn,x)
        a.len == 1 && return DA(c0)
        q0,q1,q2,p = $(fn in (:asin,:acos) ? :((u-x*x,-2x,-u,-u/2)) :
            fn == :atan ? :((u+x*x,2x,u,-u)) : fn == :asinh ? :((u+x*x,2x,u,-u/2)) :
            fn == :acosh ? :((x*x-u,2x,u,-u/2)) : :((u-x*x,-2x,-u,-u)))
        derivative = quadratic_power(q0,q1,q2,p,a.context.cutoff-1)
        $(fn == :acos) && (derivative .*= -one(eltype(derivative)))
        integral_series(a,c0,derivative)
    end
end
function Base.atan(y::DA,x::DA)
    compatible(y,x)
    x0,y0 = cons(x),cons(y)
    iszero(x0) && iszero(y0) && throw(DAError("atan is not analytic at the origin"))
    p = abs(x0) >= abs(y0) ? atan(y/x) : -atan(x/y)
    p + (atan(y0,x0)-cons(p))
end
Base.atan(a::DA,b::Real) = atan(promote(a,b)...)
Base.atan(a::Real,b::DA) = atan(promote(a,b)...)
function Base.hypot(a::DA,b::DA)
    compatible(a,b)
    scale = max(abs(cons(a)),abs(cons(b)))
    iszero(scale) && return sqrt(a*a+b*b)
    # Normalize before squaring to avoid overflow/underflow at finite centers.
    scale*sqrt((a/scale)^2+(b/scale)^2)
end
Base.hypot(a::DA,b::Real) = hypot(promote(a,b)...)
Base.hypot(a::Real,b::DA) = hypot(promote(a,b)...)
for fn in (:round,:trunc)
    @eval function Base.$fn(a::DA)
        out = copy(a)
        out.coeffs[1] = $fn(cons(a))
        finish!(out,out.len)
    end
end
function Base.mod(a::DA,p::Real)
    out = DA{promote_type(coefftype(a),typeof(p))}(a)
    out.coeffs[1] = mod(cons(a),p)
    finish!(out,out.len)
end

for fn in (:erf,:erfc)
    @eval const $fn = SpecialFunctions.$fn
    @eval function SpecialFunctions.$fn(a::DA)
        x = cons(a); c0 = SpecialFunctions.$fn(x)
        n = a.context.cutoff
        derivative = zeros(typeof(c0),n)
        derivative[1] = $(fn == :erf ? 2 : -2)*exp(-x*x)/sqrt(convert(typeof(c0),π))
        @inbounds for k in 1:n-1
            derivative[k+1] = (-2x*derivative[k] - (k > 1 ? 2derivative[k-1] : zero(c0)))/k
        end
        integral_series(a,c0,derivative)
    end
end
const loggamma = SpecialFunctions.loggamma
const gamma = SpecialFunctions.gamma
function derivative_series(f,a::DA,c0)
    a.len == 1 && return DA(c0)
    x = cons(a)
    c = Vector{typeof(c0)}(undef,a.context.cutoff+1); c[1] = c0
    factor = one(c0)
    for k in 1:length(c)-1
        factor /= k
        c[k+1] = f(k,x)*factor
    end
    series(a,c)
end
gamma_series(a::DA,c0) = derivative_series((k,x)->scalar_psi(k-1,x),a,c0)
SpecialFunctions.loggamma(a::DA) = gamma_series(a,domain_call(SpecialFunctions.loggamma,cons(a)))
function SpecialFunctions.gamma(a::DA)
    value,sign = domain_call(SpecialFunctions.logabsgamma,cons(a))
    sign*exp(gamma_series(a,value))
end
function PsiFunction(a::DA,n::Integer)
    n >= 0 || throw(ArgumentError("Invalid polygamma order"))
    x = cons(a); c0 = scalar_psi(n,x)
    derivative_series((k,x)->scalar_psi(n+k,x),a,c0)
end
SpecialFunctions.polygamma(n::Integer,a::DA) = PsiFunction(a,n)
SpecialFunctions.digamma(a::DA) = PsiFunction(a,0)

# Bessel's second-order differential equation supplies all higher derivatives.
# With y=exp(s*x)*z, scaled I/K remain scaled throughout (no overflow-prone unscaling).
function bessel_series(fn,n::Integer,a::DA,q::Int,s)
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Bessel order out of range"))
    x = cons(a); c0 = scalar_bessel(fn,n,x)
    a.len == 1 && return DA(c0)
    N = a.context.cutoff
    if iszero(x)
        fn in (SpecialFunctions.besselj,SpecialFunctions.besseli) || throw(DAError("Bessel expansion is singular at zero"))
        c = zeros(typeof(c0),N+1)
        order = abs(n)
        if order <= N
            sign = n < 0 && fn === SpecialFunctions.besselj && isodd(order) ? -1 : 1
            value = convert(typeof(c0),sign)/convert(typeof(c0),big(2)^order*factorial(big(order)))
            for k in order:2:N
                c[k+1] = value
                m = (k-order)÷2+1
                value *= -q/(convert(typeof(c0),4)*m*(m+order))
            end
        end
        return series(a,c)
    end
    c = zeros(typeof(c0),N+1); c[1] = c0
    lower,upper = scalar_bessel(fn,n-1,x),scalar_bessel(fn,n+1,x)
    neighbor = fn in (SpecialFunctions.besselj,SpecialFunctions.bessely) ? (lower-upper)/2 :
        fn in (SpecialFunctions.besseli,SpecialFunctions.besselix) ? (lower+upper)/2 : -(lower+upper)/2
    c[2] = neighbor-s*c0
    r = s*s+q
    @inbounds for k in 0:N-2
        value = -(2x*k+x+2s*x*x)*(k+1)*c[k+2] - (k*k+4s*x*k+r*x*x+s*x-n*n)*c[k+1]
        k >= 1 && (value -= (2s*(k-1)+2r*x+s)*c[k])
        k >= 2 && (value -= r*c[k-1])
        c[k+3] = value/(x*x*(k+2)*(k+1))
    end
    series(a,c)
end
for (fn,q) in ((:besselj,1),(:bessely,1),(:besseli,-1),(:besselk,-1),(:besselix,-1),(:besselkx,-1))
    @eval const $fn = SpecialFunctions.$fn
    @eval SpecialFunctions.$fn(n::Integer,a::DA) = bessel_series(SpecialFunctions.$fn,n,a,$q,
        $(fn == :besselix ? :(sign(cons(a))) : fn == :besselkx ? :(-one(cons(a))) : :(zero(cons(a)))))
end
