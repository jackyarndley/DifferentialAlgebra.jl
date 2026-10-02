"""
    adaptive_flow(advance, initial, lower, upper, times; kwargs...)

Build a [`PiecewiseTaylorMap`](@ref) of a flow at `last(times)`, splitting the
initial parameter domain when an error indicator fails at a checkpoint.
`initial(parameters)` returns the initial scalar or vector state;
`advance(state, (tstart, tend))` propagates it over one interval and returns
the same shape. Both callbacks must work with numbers and Taylor polynomials.
`times` is a finite, strictly monotone vector, with at least two entries;
backward propagation is supported. An adaptive ODE solve can implement `advance`.

All keywords and estimators match [`adaptive_map`](@ref). The initial map and
every checkpoint are checked, including independent numeric trajectories when
`check_points = true`. Guard coefficients are retained between intervals.
When a check fails, the children restart from `initial` at `first(times)`;
this recovers terms lost by the parent's truncated propagation. Accepted
intervals on an unsplit patch reuse the current state.

The returned map contains **final-time** states, even for unresolved patches
with `strict = false`. A patch is converged only if all checkpoints passed;
its `error_estimate` records the maximum absolute indicator over checkpoints.
Checks occur at the supplied times, not at the ODE solver's internal steps.
Choose checkpoint spacing and solver tolerances independently of ADS tolerances.
These heuristic checks do not certify continuous-time or uniform spatial error.
"""
function adaptive_flow(advance, initial, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}, times::AbstractVector{<:Real}; kwargs...)
    Base.require_one_based_indexing(times)
    length(times) >= 2 && all(isfinite, times) || throw(ArgumentError("Provide at least two finite checkpoint times"))
    (
        all(i -> times[i] < times[i + 1], 1:(length(times) - 1)) ||
            all(i -> times[i] > times[i + 1], 1:(length(times) - 1))
    ) ||
        throw(ArgumentError("Checkpoint times must be strictly monotone"))
    return ads_construct(CheckpointFlow(advance, initial, collect(times)), lower, upper; kwargs...)
end

struct CheckpointFlow{F, I, V}
    advance::F
    initial::I
    times::V
end

"""
    adaptive_flow(advance, initial, lower, upper, tspan::Tuple; kwargs...)

Build a flow map with online ADS checks through a solver callback.
`advance(state, (tstart, tend), monitor)` must call `monitor(state, time)`
after accepted integration steps, stop when it returns `true`, and return a
named tuple `(; state, time)` at the stopping time. `false` means continue.
The returned state must be the state passed to the last monitor call, or
the final state at `tend`. The initial map is checked before integration.

Use this overload for accepted-step monitoring; pass a vector of times to
use checkpoint monitoring with a two-argument `advance`. All other options
and the restart-on-split policy are the same. With point checks enabled,
independent numeric trajectories are advanced to each monitored time; this
can cost substantially more than coefficient-only monitoring. Validate the
final map on independent trajectories whichever estimator is chosen.
"""
function adaptive_flow(advance, initial, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}, tspan::Tuple{Real, Real}; kwargs...)
    all(isfinite, tspan) && tspan[1] != tspan[2] || throw(ArgumentError("Provide distinct finite start and end times"))
    return ads_construct(MonitoredFlow(advance, initial, promote(float.(tspan)...)), lower, upper; kwargs...)
end

struct MonitoredFlow{F, I, T}
    advance::F
    initial::I
    tspan::Tuple{T, T}
end

ads_continue(state, time) = false

function ads_flow_result(result, time)
    result isa NamedTuple && haskey(result, :state) && haskey(result, :time) ||
        throw(ArgumentError("An online propagator must return (; state, time)"))
    result.time == time || throw(ArgumentError("The propagator returned an unexpected stopping time"))
    return result.state
end

function ads_same_shape(p, first_patch)
    p.scalar == first_patch.scalar && length(p.errors) == length(first_patch.errors) ||
        throw(DimensionMismatch("The flow changed output shape"))
    eltype(p.errors) === eltype(first_patch.errors) || throw(ArgumentError("The flow changed coefficient type"))
    return nothing
end

function ads_flow_check(value, box, probes, states, ctx, options)
    p = ads_analyze(value, box, options, ctx)
    options.check_points && ads_check_points!(p, probes, states)
    ads_check_context(ctx)
    return ads_assess(p, options)
end

function ads_candidate(problem::CheckpointFlow, box, x, ctx, options, allow_split)
    state = problem.initial(box.center .+ box.radius .* x)
    probes = options.check_points ? ads_probes(box) : nothing
    numeric_states = options.check_points ? [problem.initial(point) for (point, _, _) in probes] : nothing
    first_patch = ads_flow_check(state, box, probes, numeric_states, ctx, options)
    p = first_patch
    errors = copy(p.errors)
    accepted = p.accepted
    for k in 1:(length(problem.times) - 1)
        # Never translate a failed, already truncated trajectory into children.
        !p.accepted && allow_split && return p
        span = (problem.times[k], problem.times[k + 1])
        state = problem.advance(state, span)
        if options.check_points
            numeric_states = [problem.advance(u, span) for u in numeric_states]
        end
        p = ads_flow_check(state, box, probes, numeric_states, ctx, options)
        ads_same_shape(p, first_patch)
        errors .= max.(errors, p.errors)
        accepted &= p.accepted
    end
    return (; p..., errors, accepted)
end

function ads_candidate(problem::MonitoredFlow, box, x, ctx, options, allow_split)
    state = problem.initial(box.center .+ box.radius .* x)
    probes = options.check_points ? ads_probes(box) : nothing
    numeric_states = options.check_points ? [problem.initial(point) for (point, _, _) in probes] : nothing
    first_patch = ads_flow_check(state, box, probes, numeric_states, ctx, options)
    !first_patch.accepted && allow_split && return first_patch
    errors = copy(first_patch.errors)
    accepted = Ref(first_patch.accepted)
    current = Ref(first_patch)
    previous_time = Ref(first(problem.tspan))
    stopped = Ref(false)
    direction = sign(last(problem.tspan) - first(problem.tspan))

    function monitor(value, time)
        stopped[] && throw(ArgumentError("The propagator continued after the monitor requested a stop"))
        isfinite(time) && direction * (time - previous_time[]) >= 0 &&
            direction * (last(problem.tspan) - time) >= 0 ||
            throw(ArgumentError("Monitor times must progress toward the end of tspan"))
        if options.check_points && time != previous_time[]
            span = (previous_time[], time)
            for i in eachindex(numeric_states)
                numeric_states[i] = ads_flow_result(problem.advance(numeric_states[i], span, ads_continue), time)
            end
        end
        p = ads_flow_check(value, box, probes, numeric_states, ctx, options)
        ads_same_shape(p, first_patch)
        errors .= max.(errors, p.errors)
        accepted[] &= p.accepted
        current[] = p
        previous_time[] = time
        stopped[] = !p.accepted && allow_split
        return stopped[]
    end

    result = problem.advance(state, problem.tspan, monitor)
    target = stopped[] ? previous_time[] : last(problem.tspan)
    final = ads_flow_result(result, target)
    # A solver need not invoke its callback at the final endpoint.
    if !stopped[]
        monitor(final, target)
    end
    p = current[]
    return (; p..., errors, accepted = accepted[])
end
