# Julia's generic normalize uses floating-point limits to choose a rescaling
# strategy. A polynomial has no typemax/prevfloat; divide by its analytic norm.
function LinearAlgebra.normalize(a::AbstractArray{<:TaylorPolynomial}, p::Real = 2)
    return isempty(a) ? copy(a) : a ./ LinearAlgebra.norm(a, p)
end
function LinearAlgebra.normalize!(a::AbstractArray{<:TaylorPolynomial}, p::Real = 2)
    isempty(a) && return a
    a ./= LinearAlgebra.norm(a, p)
    return a
end

scalar_eigen_type(::Type{T}) where {T} =
    promote_type(Float32, typeof(zero(T) / sqrt(abs2(one(T)))))

"""
    eigenpairs(A)

Eigenvalues and orthonormal eigenvectors of a real symmetric polynomial matrix.
Diagonal matrices use coordinate eigenvectors, sorted by constant eigenvalues
with stable ties. Other nonconstant matrices require distinct constant eigenvalues
for a unique Taylor expansion.
Each eigenpair is lifted degree by degree with the constant bordered Jacobian.
"""
function eigenpairs(A::AbstractMatrix{<:TaylorPolynomial})
    Base.require_one_based_indexing(A)
    n, m = size(A)
    n == m || throw(DimensionMismatch("Expected a square matrix"))
    if n == 0
        T = isconcretetype(eltype(A)) ? coefficient_type(eltype(A)) : Float64
        # Match the scalar solver's real working type without an empty solve.
        R = scalar_eigen_type(T)
        return TaylorPolynomial{R}[], Matrix{TaylorPolynomial{R}}(undef, 0, 0)
    end
    ctx = compatible(A)
    limit = ctx.basis.ends[working_order(ctx) + 1]
    for j in 1:n, i in 1:(j - 1)
        a, b = A[i, j], A[j, i]
        stop = min(limit, max(a.len, b.len))
        all(k -> coefficient_at(a, k) == coefficient_at(b, k), 1:stop) ||
            throw(ArgumentError("Expected a symmetric matrix"))
    end
    T = mapreduce(coefficient_type, promote_type, A)
    if all(
            i == j || all(
                coefficient_iszero,
                @view(A[i, j].coeffs[1:min(A[i, j].len, limit)])
            ) for i in 1:n, j in 1:n
        )
        R = scalar_eigen_type(T)
        diagonal = LinearAlgebra.diag(A)
        permutation = sortperm(constant_term.(diagonal))
        values = TaylorPolynomial{R}.(diagonal[permutation])
        for value in values
            finish!(value, min(value.len, limit))
        end
        vectors = TaylorPolynomial{R}.(Ref(ctx), Matrix{R}(I, n, n)[:, permutation])
        return values, vectors
    end
    A0 = Matrix{T}(constant_term.(A))
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(A0))
    R, V = eltype(F.values), eltype(F.vectors)
    values = TaylorPolynomial{R}.(Ref(ctx), F.values; filter = false)
    vectors = TaylorPolynomial{V}.(Ref(ctx), F.vectors; filter = false)
    all(a -> all(coefficient_iszero, @view(a.coeffs[2:min(a.len, limit)])), A) &&
        return finish!.(values), finish!.(vectors)
    scale = max(one(R), maximum(abs, F.values; init = zero(R)))
    minimum((abs(F.values[i] - F.values[j]) for i in 1:n for j in (i + 1):n); init = convert(R, Inf)) > 100eps(R) * scale ||
        throw(ArgumentError("Repeated or numerically indistinguishable eigenvalues do not define unique Taylor eigenvectors"))
    for k in 1:n
        v0, lambda0 = F.vectors[:, k], F.values[k]
        J = [A0 - lambda0 * I -v0; permutedims(v0) zero(R)]
        Ji = inv(J)
        v, lambda = copy.(vectors[:, k]), copy(values[k])
        for degree in 1:working_order(ctx)
            residual = eigen_residual(A, v, lambda, degree)
            correction = _linear_transform(Ji, residual; filter = false)
            for i in eachindex(v)
                weighted_sum!(v[i], v[i], one(V), correction[i], -one(V); filter = false)
            end
            weighted_sum!(lambda, lambda, one(R), correction[end], -one(R); filter = false)
        end
        vectors[:, k] = v
        values[k] = lambda
    end
    return finish!.(values), finish!.(vectors)
end

# Local degree limits make eigenpair lifting safe alongside ordinary arithmetic.
# Reuse one product buffer instead of allocating every scalar matrix product.
function eigen_residual(A, v, lambda, degree::Int)
    ctx = valid(lambda); T = coefficient_type(lambda)
    n = length(v)
    result = [allocate(ctx, T, ctx.basis.ends[degree + 1]) for _ in 1:(n + 1)]
    work = allocate(ctx, T, ctx.basis.ends[degree + 1])
    for i in 1:n
        for j in 1:n
            multiply!(work, A[i, j], v[j], degree; filter = false)
            weighted_sum!(result[i], result[i], one(T), work, one(T); order = degree, filter = false)
        end
        multiply!(work, lambda, v[i], degree; filter = false)
        weighted_sum!(result[i], result[i], one(T), work, -one(T); order = degree, filter = false)
        multiply!(work, v[i], v[i], degree; filter = false)
        weighted_sum!(result[end], result[end], one(T), work, one(T) / 2; order = degree, filter = false)
    end
    result[end].coeffs[1] -= one(T) / 2
    finish!(result[end], result[end].len; filter = false)
    return result
end
function eigenpairs(A::AbstractMatrix{<:Real})
    Base.require_one_based_indexing(A)
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("Expected a square matrix"))
    LinearAlgebra.issymmetric(A) || throw(ArgumentError("Expected a symmetric matrix"))
    F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(A)))
    return F.values, F.vectors
end
