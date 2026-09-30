"All multi-indices of total degree <= order, ordered by degree then reverse lexicographically."
function getMultiIndices(no::Integer, nv::Integer)
    0 <= no <= 65535 && 1 <= nv <= 1024 || throw(ArgumentError("Invalid order or dimension"))
    count = binomial(big(no) + nv, nv)
    count * nv <= 32 * 1024^2 || throw(ArgumentError("Monomial basis is too large"))
    result = Vector{UInt32}[]
    sizehint!(result, Int(count))
    foreach_monomial(Int(no), Int(nv)) do powers, _
        push!(result, UInt32.(powers))
    end
    return result
end
function getRawMoments(mgf::DA, no::Integer)
    return lock(engine_lock) do
        valid(mgf)
        0 <= no <= getMaxOrder() || throw(ArgumentError("Moment order exceeds initialized order"))
        indices = getMultiIndices(no, getMaxVariables())
        T = coefftype(mgf)
        moments = [getCoefficient(mgf, jj) * prod(j -> convert(T, factorial(big(j))), jj) for jj in indices]
        indices, moments
    end
end
function getCentralMoments(mgf::DA, no::Integer)
    return lock(engine_lock) do
        valid(mgf)
        mean = linear(mgf)
        shift = sum(mean[i] * variable(i, coefftype(mgf)) for i in eachindex(mean))
        getRawMoments(exp(-shift) * mgf, no)
    end
end
