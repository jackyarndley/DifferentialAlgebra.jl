# Named right-hand sides share solver specializations across regression cases.
oscillator_state(u, p, t) = [u[2], -u[1]]
function oscillator_state!(du, u, p, t)
    du[1], du[2] = u[2], -u[1]
    return nothing
end
exponential_growth(u, p, t) = u
function exponential_growth!(du, u, p, t)
    du[1] = u[1]
    return nothing
end

function zeroed_oscillator!(du, u, strategy, t)
    if strategy === :fill
        fill!(du, 0.0)
    else
        du .= 0
    end
    du[1] += u[2]
    du[2] -= u[1]
    return nothing
end

function recorded_oscillator(u, strategy, t)
    du = DifferentialAlgebra.time_rhs_buffer(u)
    zeroed_oscillator!(du, u, strategy, t)
    return du
end

parameter_rate(p) = p isa Real ? p : first(p)
parameter_growth(u, p, t) = parameter_rate(p) .* u
function parameter_growth!(du, u, p, t)
    du[1] = parameter_rate(p) * u[1]
    return nothing
end
