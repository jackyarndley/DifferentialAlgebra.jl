# Named right-hand sides share solver specializations across regression cases.
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
