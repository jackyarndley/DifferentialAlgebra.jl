# Closed scalar derivatives at higher precision, independent of the model's
# normalized composition and interval remainder formulas.
function scalar_model_series(f, center, order)
    return setprecision(768) do
        c = BigFloat(center)
        map(0:order) do k
            k == 0 && return f(c)
            f === inv && return (-1)^k / c^(k + 1)
            f === log && return (-1)^(k + 1) / (k * c^k)
            factor = prod(big(1) // 2 - j for j in 0:(k - 1)) / factorial(big(k))
            factor * sqrt(c) / c^k
        end
    end
end

function scalar_model_residual(f, series, center, displacement)
    return setprecision(768) do
        c, h = BigFloat(center), BigFloat(displacement)
        approximation = sum(series[k + 1] * h^k for k in 0:(length(series) - 1))
        value = f(c + h)
        value, value - approximation
    end
end
