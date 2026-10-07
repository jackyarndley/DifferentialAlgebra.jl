module IntervalTestSupport
using DifferentialAlgebra, IntervalArithmetic

const IA = IntervalArithmetic
const DA = DifferentialAlgebra
const IF = Interval{Float64}
const IB = Interval{BigFloat}
subset(a, b) = IA.issubset_interval(a, b)
interval_contains(a, x) = IA.in_interval(x, a)
sameinterval(a, b) = IA.isequal_interval(a, b)
guaranteed_coefficients(p) = all(IA.isguaranteed, p.coeffs[1:p.len])
end
