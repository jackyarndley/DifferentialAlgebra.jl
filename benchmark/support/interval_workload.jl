# Keep the six-variable expression identical across interval benchmarks so that
# construction, arithmetic, enclosure and ADS measurements stay comparable.
uncertainty_map(x) = exp(x[1] + x[2] / 2) * cos(x[3] - x[4]) +
    log(2 + x[5] * x[6])

bound_summary(x) = (;
    lower = IntervalArithmetic.inf(x), upper = IntervalArithmetic.sup(x),
    decoration = IntervalArithmetic.decoration(x),
    guaranteed = IntervalArithmetic.isguaranteed(x),
)
