# Exact allocation counts apply to supported, uninstrumented stable runtimes.
# Numerical, aliasing and ownership checks run on coverage and prerelease jobs too.
const CHECK_ZERO_ALLOCATIONS = iszero(Base.JLOptions().code_coverage) &&
    isempty(VERSION.prerelease)
