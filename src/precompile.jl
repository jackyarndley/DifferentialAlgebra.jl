using PrecompileTools: @compile_workload

# One small Float64 workload covers the common scalar kernels. Dimension and
# order are runtime data, so this also helps larger polynomial bases. Keep
# integrators, plotting and user callbacks out of the
# package cache. Never serialize an initialized global algebra.
@compile_workload begin
    with_algebra(4, 2) do _
        x, y = variable(1), variable(2)
        p = exp(x / 10) * sin(y / 10) + sqrt(2 + x - y) / (3 + x + y)
        map = CompiledMap([p, differentiate(p, 1)])
        evaluate(map, [0.1, 0.2])
    end
end
