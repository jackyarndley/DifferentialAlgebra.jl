# Record an analytic RHS once, then compute one new coefficient per node and
# degree. Runtime graph storage keeps compilation independent of graph size.
@enum TimeOperation::UInt8 begin
    TimeInput
    TimeConstant
    TimeClock
    TimeAdd
    TimeSubtract
    TimeMultiply
    TimeDivide
    TimePower
    TimeExp
    TimeLog
    TimeSin
    TimeCos
    TimeSinh
    TimeCosh
end
struct TimeInstruction{R}
    operation::TimeOperation
    left::Int
    right::Int
    power::R
end
struct TimeTape{T, R}
    instructions::Vector{TimeInstruction{R}}
    constants::Vector{T}
end
struct TimeNode{T, R} <: Number
    tape::TimeTape{T, R}
    index::Int
end
function time_node(tape::TimeTape{T, R}, op, value, left = 0, right = 0, power = zero(R)) where {T, R}
    isfinite(value) || throw(DomainError(value, "Nonfinite time coefficient"))
    push!(tape.instructions, TimeInstruction{R}(op, left, right, power))
    push!(tape.constants, value)
    return TimeNode(tape, length(tape.constants))
end
time_value(a::TimeNode) = a.tape.constants[a.index]
time_constant(a::TimeNode, x::Real) = time_node(a.tape, TimeConstant, x)
Base.zero(a::TimeNode) = time_constant(a, zero(time_value(a)))
Base.one(a::TimeNode) = time_constant(a, one(time_value(a)))
Base.copy(a::TimeNode) = a
Base.real(a::TimeNode) = a
Base.conj(a::TimeNode) = a
Base.abs2(a::TimeNode) = a * a
Base.:+(a::TimeNode) = a
Base.:-(a::TimeNode) = zero(a) - a
function record_time_binary(op, f, a::TimeNode, b::TimeNode)
    a.tape === b.tape || throw(ArgumentError("Cannot mix independent time expansions"))
    return time_node(a.tape, op, f(time_value(a), time_value(b)), a.index, b.index)
end
for (f, op) in ((:+, :TimeAdd), (:-, :TimeSubtract), (:*, :TimeMultiply), (:/, :TimeDivide))
    @eval begin
        Base.$f(a::TimeNode, b::TimeNode) = record_time_binary($op, $f, a, b)
        Base.$f(a::TimeNode, b::Real) = $f(a, time_constant(a, b))
        Base.$f(a::Real, b::TimeNode) = $f(time_constant(b, a), b)
    end
end
Base.inv(a::TimeNode) = one(a) / a
function Base.:^(a::TimeNode, n::Integer)
    n < 0 && return inv(a)^(-(n + 1)) / a
    n == 0 && return one(a)
    n == 1 && return a
    half = a^(n >> 1)
    square = half * half
    return isodd(n) ? square * a : square
end
function record_time_power(a::TimeNode, p::Real)
    p isa TaylorPolynomial && return exp(p * log(a))
    isinteger(p) && return a^BigInt(p)
    return time_node(a.tape, TimePower, time_value(a)^p, a.index, 0, p)
end
Base.:^(a::TimeNode, p::Real) = record_time_power(a, p)
Base.:^(a::TimeNode, p::Rational) = record_time_power(a, p)
Base.:^(a::TimeNode, b::TimeNode) = exp(b * log(a))
Base.:^(a::Real, b::TimeNode) = exp(log(a) * b)
Base.:^(::Irrational{:ℯ}, b::TimeNode) = exp(b)
Base.sqrt(a::TimeNode) = time_node(a.tape, TimePower, sqrt(time_value(a)), a.index, 0, 1 // 2)
Base.cbrt(a::TimeNode) = time_node(a.tape, TimePower, cbrt(time_value(a)), a.index, 0, 1 // 3)
Base.exp(a::TimeNode) = time_node(a.tape, TimeExp, exp(time_value(a)), a.index)
Base.log(a::TimeNode) = time_node(a.tape, TimeLog, log(time_value(a)), a.index)
function record_time_sincos(a::TimeNode, hyperbolic)
    x = time_value(a)
    s, c = hyperbolic ? (sinh(x), cosh(x)) : sincos(x)
    i = length(a.tape.constants) + 1
    snode = time_node(a.tape, hyperbolic ? TimeSinh : TimeSin, s, a.index, i + 1)
    cnode = time_node(a.tape, hyperbolic ? TimeCosh : TimeCos, c, a.index, i)
    return snode, cnode
end
Base.sincos(a::TimeNode) = record_time_sincos(a, false)
for (f, hyperbolic, index) in ((:sin, false, 1), (:cos, false, 2), (:sinh, true, 1), (:cosh, true, 2))
    @eval Base.$f(a::TimeNode) = record_time_sincos(a, $hyperbolic)[$index]
end
Base.tan(a::TimeNode) = sin(a) / cos(a)
Base.tanh(a::TimeNode) = sinh(a) / cosh(a)

function time_convolution(a, b, k, prototype)
    value = zero(prototype)
    @inbounds for j in 0:k
        aj, bj = a[j + 1], b[k - j + 1]
        (iszero(aj) || iszero(bj)) && continue
        value = muladd(aj, bj, value)
    end
    return value
end
function time_convolution(a::Vector{T}, b::Vector{T}, k, prototype::T) where {T <: Union{Float32, Float64}}
    value = zero(T)
    @inbounds @simd for j in 0:k
        value += a[j + 1] * b[k - j + 1]
    end
    return value
end
function time_convolution(a::Vector{P}, b::Vector{P}, k, prototype::P) where {P <: TaylorPolynomial}
    value, scratch = zero(prototype), zero(prototype)
    @inbounds for j in 0:k
        (iszero(a[j + 1]) || iszero(b[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, a[j + 1], b[k - j + 1])
        add!(value, value, scratch)
    end
    return value
end

@inline time_weight(j, k, p, ::Val{:exp}) = j
@inline time_weight(j, k, p, ::Val{:log}) = j - k
@inline time_weight(j, k, p, ::Val{:power}) = p * j - (k - j)
function time_weighted_sum(a, c, k, p, mode)
    value = zero(c[1])
    @inbounds for j in 1:k
        weight = time_weight(j, k, p, mode)
        iszero(weight) && continue
        value += weight * a[j + 1] * c[k - j + 1]
    end
    return value
end
function time_weighted_sum(a::Vector{T}, c::Vector{T}, k, p, mode) where {T <: Union{Float32, Float64}}
    value = zero(T)
    @inbounds @simd for j in 1:k
        value += time_weight(j, k, p, mode) * a[j + 1] * c[k - j + 1]
    end
    return value
end
function time_weighted_sum(a::Vector{P}, c::Vector{P}, k, p, mode) where {P <: TaylorPolynomial}
    value, scratch = zero(c[1]), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        weight = time_weight(j, k, p, mode)
        (iszero(weight) || iszero(a[j + 1]) || iszero(c[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, a[j + 1], c[k - j + 1])
        weighted_sum!(value, value, unit, scratch, weight)
    end
    return value
end
time_normalize(value, reciprocal, k) = (value / k) * reciprocal
function time_normalize(value::TaylorPolynomial, reciprocal::TaylorPolynomial, k)
    scale!(value, value, inv(convert(coefficient_type(value), k)))
    return value * reciprocal
end
function time_quotient_coefficient(a, b, c, k, reciprocal)
    value = a[k + 1]
    @inbounds for j in 1:k
        value -= b[j + 1] * c[k - j + 1]
    end
    return time_normalize(value, reciprocal, 1)
end
function time_quotient_coefficient(a::Vector{P}, b::Vector{P}, c::Vector{P}, k, reciprocal) where {P <: TaylorPolynomial}
    value, scratch = copy(a[k + 1]), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        (iszero(b[j + 1]) || iszero(c[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, b[j + 1], c[k - j + 1])
        weighted_sum!(value, value, unit, scratch, -unit)
    end
    return time_normalize(value, reciprocal, 1)
end

function advance_time_coefficient!(values, instructions, reciprocals, k)
    @inbounds for node in eachindex(instructions)
        instruction = instructions[node]
        op, left, right = instruction.operation, instruction.left, instruction.right
        c = values[node]
        op === TimeInput && continue # Updated from the ODE after each degree.
        if op === TimeConstant
            c[k + 1] = zero(c[1])
        elseif op === TimeClock
            c[k + 1] = k == 1 ? one(c[1]) : zero(c[1])
        elseif op === TimeAdd
            c[k + 1] = values[left][k + 1] + values[right][k + 1]
        elseif op === TimeSubtract
            c[k + 1] = values[left][k + 1] - values[right][k + 1]
        elseif op === TimeMultiply
            a, b = values[left], values[right]
            if instructions[left].operation === TimeConstant
                c[k + 1] = a[1] * b[k + 1]
            elseif instructions[right].operation === TimeConstant
                c[k + 1] = a[k + 1] * b[1]
            else
                c[k + 1] = time_convolution(a, b, k, c[1])
            end
        elseif op === TimeDivide
            a, b = values[left], values[right]
            c[k + 1] = instructions[right].operation === TimeConstant ? a[k + 1] * reciprocals[node] : time_quotient_coefficient(a, b, c, k, reciprocals[node])
        elseif op === TimeExp
            a = values[left]
            c[k + 1] = time_weighted_sum(a, c, k, 0, Val(:exp)) / k
        elseif op === TimePower || op === TimeLog
            a = values[left]
            iszero(constant_term(a[1])) && throw(DomainError(a[1], "The time expansion is not analytic here"))
            p = instruction.power
            value = op === TimeLog ? time_weighted_sum(a, c, k, p, Val(:log)) : time_weighted_sum(a, c, k, p, Val(:power))
            op === TimeLog && (value += k * a[k + 1])
            c[k + 1] = time_normalize(value, reciprocals[node], k)
        else # Coupled trigonometric recurrences use only earlier degrees.
            a, partner = values[left], values[right]
            value = time_weighted_sum(a, partner, k, 0, Val(:exp))
            c[k + 1] = (op === TimeCos ? -value : value) / k
        end
    end
    return nothing
end

struct TimeWorkspace{T, R}
    tape::TimeTape{T, R}
    values::Vector{Vector{T}}
    states::Vector{TimeNode{T, R}}
    outputs::Vector{Int}
    series::Vector{TimeSeries{T}}
    reciprocals::Vector{T}
    order::Int
    nstate::Int
end
function time_workspace(f, initial, epoch, order, parameters, ::Type{T}) where {T <: Real}
    R = typeof(float(constant_term(zero(T))))
    tape = TimeTape(TimeInstruction{R}[], T[])
    workspace = TimeWorkspace(tape, Vector{T}[], TimeNode{T, R}[], Int[], TimeSeries{T}[], T[], order, length(initial))
    update_time_expansion!(workspace, f, initial, epoch, parameters)
    return workspace
end

# Re-record at every step, respecting parameter/callback changes while reusing
# node and coefficient buffers. No generated method or fixed graph is retained.
function update_time_expansion!(workspace::TimeWorkspace{T}, f, initial, epoch, parameters) where {T}
    (; tape, values, states, outputs, series, reciprocals, order, nstate) = workspace
    length(initial) == nstate || throw(DimensionMismatch("The state dimension changed"))
    empty!(tape.instructions)
    empty!(tape.constants)
    empty!(states)
    empty!(outputs)
    for value in initial
        push!(states, time_node(tape, TimeInput, value))
    end
    time = time_node(tape, TimeClock, epoch)
    derivative = f(states, parameters, time)
    derivative isa AbstractVector && length(derivative) == nstate ||
        throw(DimensionMismatch("The right-hand side must return one derivative per state"))
    for value in derivative
        if value isa TimeNode
            value.tape === tape || throw(ArgumentError("Cannot mix independent time expansions"))
            push!(outputs, value.index)
        else
            value isa Real || throw(ArgumentError("Unsupported right-hand-side element"))
            push!(outputs, time_node(tape, TimeConstant, value).index)
        end
    end
    while length(values) < length(tape.constants)
        push!(values, Vector{T}(undef, order + 1))
    end
    for i in eachindex(tape.constants)
        values[i][1] = copy(tape.constants[i])
    end
    resize!(reciprocals, length(tape.instructions))
    for i in eachindex(tape.instructions)
        node = tape.instructions[i]
        if node.operation in (TimeDivide, TimePower, TimeLog)
            denominator = values[node.operation === TimeDivide ? node.right : node.left][1]
            iszero(constant_term(denominator)) && throw(DomainError(denominator, "The time expansion is not analytic here"))
            reciprocals[i] = inv(denominator)
        end
    end
    for k in 0:(order - 1)
        k > 0 && advance_time_coefficient!(values, tape.instructions, reciprocals, k)
        for i in 1:nstate
            value = values[outputs[i]][k + 1] / (k + 1)
            isfinite(value) || throw(DomainError(value, "Nonfinite time coefficient"))
            values[i][k + 2] = value
        end
    end
    if isempty(series)
        for i in 1:nstate
            push!(series, TimeSeries{T}(values[i]))
        end
    end
    return series
end
