# Record an analytic RHS once per expansion, then compute one new coefficient
# per node and degree using time_kernels.jl. Runtime graph storage keeps
# compilation independent of graph size.
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
# An in-place RHS can assign ordinary constants while recording a time graph.
# Keep the storage concrete and attach those constants to the same tape.
struct TimeDerivativeBuffer{N <: TimeNode} <: AbstractVector{N}
    data::Vector{N}
end
Base.size(buffer::TimeDerivativeBuffer) = size(buffer.data)
Base.IndexStyle(::Type{<:TimeDerivativeBuffer}) = IndexLinear()
Base.getindex(buffer::TimeDerivativeBuffer, i::Int) = buffer.data[i]
function Base.setindex!(buffer::TimeDerivativeBuffer, value, i::Int)
    checkbounds(buffer.data, i)
    reference = first(buffer.data)
    node = value isa TimeNode ? value : time_constant(reference, value)
    node.tape === reference.tape || throw(
        ArgumentError(
            "Cannot mix independent time expansions"
        )
    )
    buffer.data[i] = node
    return value
end
function Base.fill!(buffer::TimeDerivativeBuffer, value)
    for i in eachindex(buffer)
        buffer[i] = value
    end
    return buffer
end
time_rhs_buffer(u::AbstractVector) = zero.(u)
time_rhs_buffer(u::AbstractVector{<:TimeNode}) = TimeDerivativeBuffer(zero.(u))

function time_node(tape::TimeTape{T, R}, op, value, left = 0, right = 0, power = zero(R)) where {T, R}
    isfinite(value) || throw(TaylorError("Nonfinite time coefficient"))
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
    return time_node(
        a.tape, op, domain_call(x -> f(x, time_value(b)), time_value(a)), a.index, b.index
    )
end
for (f, op) in ((:+, :TimeAdd), (:-, :TimeSubtract), (:*, :TimeMultiply), (:/, :TimeDivide))
    @eval begin
        Base.$f(a::TimeNode, b::TimeNode) = record_time_binary($op, $f, a, b)
        Base.$f(a::TimeNode, b::Real) = $f(a, time_constant(a, b))
        Base.$f(a::Real, b::TimeNode) = $f(time_constant(b, a), b)
    end
end
for (A, B, C) in (
        (TimeNode, TimeNode, TimeNode), (TimeNode, TimeNode, Real),
        (TimeNode, Real, TimeNode), (Real, TimeNode, TimeNode),
        (TimeNode, Real, Real), (Real, TimeNode, Real), (Real, Real, TimeNode),
    )
    # Record the existing operations without promoting a node away from its tape.
    @eval Base.muladd(a::$A, b::$B, c::$C) = a * b + c
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
    return time_node(a.tape, TimePower, domain_call(x -> x^p, time_value(a)), a.index, 0, p)
end
Base.:^(a::TimeNode, p::Real) = record_time_power(a, p)
Base.:^(a::TimeNode, p::Rational) = record_time_power(a, p)
Base.:^(a::TimeNode, b::TimeNode) = exp(b * log(a))
function Base.:^(a::Real, b::TimeNode)
    base = convert(promote_type(typeof(a), typeof(time_value(b))), a)
    coefficient_iszero(base) && return time_constant(
        b, zero_base_power(
            base, time_value(b),
            b.tape.instructions[b.index].operation === TimeConstant && iszero(time_value(b))
        )
    )
    return exp(domain_call(log, base) * b)
end
Base.:^(::Irrational{:ℯ}, b::TimeNode) = exp(b)
Base.sqrt(a::TimeNode) = time_node(
    a.tape, TimePower, domain_call(sqrt, time_value(a)), a.index, 0, 1 // 2
)
Base.cbrt(a::TimeNode) = time_node(
    a.tape, TimePower, domain_call(cbrt, time_value(a)), a.index, 0, 1 // 3
)
Base.exp(a::TimeNode) = time_node(a.tape, TimeExp, exp(time_value(a)), a.index)
Base.log(a::TimeNode) = time_node(a.tape, TimeLog, domain_call(log, time_value(a)), a.index)
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
function Base.tan(a::TimeNode)
    s, c = record_time_sincos(a, false)
    return s / c
end
function Base.tanh(a::TimeNode)
    s, c = record_time_sincos(a, true)
    return s / c
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
            c[k + 1] = if instructions[right].operation === TimeConstant
                a[k + 1] * reciprocals[node]
            else
                time_quotient_coefficient(a, b, c, k, reciprocals[node])
            end
        elseif op === TimeExp
            a = values[left]
            c[k + 1] = time_unary_coefficient(a, c, k, 0, nothing, Val(:exp))
        elseif op === TimePower || op === TimeLog
            a = values[left]
            !invertible_constant(constant_term(a[1])) && throw(TaylorError("The time expansion is not analytic here"))
            mode = op === TimeLog ? Val(:log) : Val(:power)
            c[k + 1] = time_unary_coefficient(a, c, k, instruction.power, reciprocals[node], mode)
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
    workspace = TimeWorkspace(
        tape, Vector{T}[], TimeNode{T, R}[], Int[], TimeSeries{T}[], T[],
        order, length(initial)
    )
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
            input = node.operation === TimeDivide ? node.right : node.left
            denominator = values[input][1]
            !invertible_constant(constant_term(denominator)) &&
                throw(TaylorError("The time expansion is not analytic here"))
            reciprocals[i] = inv(denominator)
        end
    end
    for k in 0:(order - 1)
        k > 0 && advance_time_coefficient!(values, tape.instructions, reciprocals, k)
        for i in 1:nstate
            value = values[outputs[i]][k + 1] / (k + 1)
            isfinite(value) || throw(TaylorError("Nonfinite time coefficient"))
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
