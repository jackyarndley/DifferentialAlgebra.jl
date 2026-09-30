gradient(a::DA) = [deriv(a, i) for i in 1:getMaxVariables()]
jacobian(a::DA) = permutedims(gradient(a))
jacobian(v::AbstractVector{<:DA}) = [deriv(v[i], j) for i in eachindex(v), j in 1:getMaxVariables()]
hessian(a::DA) = [deriv(deriv(a, i), j) for i in 1:getMaxVariables(), j in 1:getMaxVariables()]
hessian(v::AbstractVector{<:DA}) = hessian.(v)

# Julia's generic normalize uses floating-point limits to choose a rescaling
# strategy. A polynomial has no typemax/prevfloat; divide by its analytic norm.
function LinearAlgebra.normalize(a::AbstractArray{<:DA}, p::Real = 2)
    return isempty(a) ? copy(a) : a ./ LinearAlgebra.norm(a, p)
end
function LinearAlgebra.normalize!(a::AbstractArray{<:DA}, p::Real = 2)
    isempty(a) && return a
    a ./= LinearAlgebra.norm(a, p)
    return a
end

function hess_stack(v::AbstractVector{<:DA})
    if isempty(v)
        T = isconcretetype(eltype(v)) ? eltype(v) : DA{Float64}
        return Array{T}(undef, getMaxVariables(), getMaxVariables(), 0)
    end
    return stack(hessian(v); dims = 3)
end

"""
    eigh(A)

Eigenvalues and orthonormal eigenvectors of a real symmetric polynomial matrix.
The constant matrix must have distinct eigenvalues for a unique Taylor expansion.
Each eigenpair is lifted degree by degree with the constant bordered Jacobian.
"""
function eigh(A::AbstractMatrix{<:DA})
    Base.require_one_based_indexing(A)
    foreach(valid, A)
    n, m = size(A)
    n == m || throw(DimensionMismatch("Expected a square matrix"))
    if n == 0
        T = isconcretetype(eltype(A)) ? coefftype(eltype(A)) : Float64
        return DA{T}[], Matrix{DA{T}}(undef, 0, 0)
    end
    T = mapreduce(coefftype, promote_type, A)
    all(norm(A[i, j] - A[j, i], 0) == 0 for i in 1:n,j in 1:n) || throw(ArgumentError("Expected a symmetric matrix"))
    if all(iszero(A[i, j]) for i in 1:n for j in 1:n if i != j)
        permutation = sortperm(cons.(LinearAlgebra.diag(A)))
        return copy.(LinearAlgebra.diag(A)[permutation]), DA{T}.(Matrix{T}(I, n, n)[:, permutation])
    end
    A0 = Matrix{T}(cons.(A))
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(A0))
    values, vectors = DA{T}.(F.values), DA{T}.(F.vectors)
    all(a -> a.len == 1, A) && return values, vectors
    scale = max(one(T), maximum(abs, F.values; init = zero(T)))
    minimum((abs(F.values[i] - F.values[j]) for i in 1:n for j in (i + 1):n); init = convert(T, Inf)) > 100eps(T) * scale ||
        throw(ArgumentError("Repeated or numerically indistinguishable eigenvalues do not define unique Taylor eigenvectors"))
    ctx = valid(first(A))
    for k in 1:n
        v0, lambda0 = F.vectors[:, k], F.values[k]
        J = [A0 - lambda0 * I -v0; permutedims(v0) zero(T)]
        Ji = inv(J)
        v, lambda = copy(vectors[:, k]), values[k]
        for degree in 1:ctx.cutoff
            residual = eigen_residual(A, v, lambda, degree)
            correction = _linear_transform(Ji, residual)
            v = v - correction[1:n]
            lambda = lambda - correction[end]
        end
        vectors[:, k] = v
        values[k] = lambda
    end
    return values, vectors
end

# Local degree limits make eigenpair lifting safe alongside ordinary arithmetic.
# Reuse one product buffer instead of allocating every scalar matrix product.
function eigen_residual(A, v, lambda, degree::Int)
    ctx = valid(lambda); T = coefftype(lambda)
    n = length(v)
    result = [allocate(ctx, T, ctx.basis.ends[degree + 1]) for _ in 1:(n + 1)]
    work = allocate(ctx, T, ctx.basis.ends[degree + 1])
    for i in 1:n
        for j in 1:n
            multiply!(work, A[i, j], v[j], degree)
            weighted_sum!(result[i], result[i], one(T), work, one(T); order = degree)
        end
        multiply!(work, lambda, v[i], degree)
        weighted_sum!(result[i], result[i], one(T), work, -one(T); order = degree)
        multiply!(work, v[i], v[i], degree)
        weighted_sum!(result[end], result[end], one(T), work, one(T) / 2; order = degree)
    end
    result[end].coeffs[1] -= one(T) / 2
    finish!(result[end], result[end].len)
    return result
end
function eigh(A::AbstractMatrix{<:Real})
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(A)))
    return F.values, F.vectors
end
