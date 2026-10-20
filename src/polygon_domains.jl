# Exact geometry is outside coefficient kernels. Stored binary floats are
# converted to exact rationals so clipping, coverage and membership need no eps.
const PolygonReal = Rational{BigInt}
const PolygonPoint = NTuple{2, PolygonReal}

function polygon_real(x::Union{Integer, Rational, AbstractFloat})
    isfinite(x) || throw(ArgumentError("Polygon coordinates must be finite"))
    return deepcopy(PolygonReal(x))
end
polygon_real(x) = throw(ArgumentError("Supply finite integer, rational or floating polygon coordinates"))

# These bounds are only a filter for exact geometry. Nonfinite arithmetic gives
# an inconclusive enclosure and therefore cannot exclude a point or support.
ads_float_bounds(values) = (
    map(x -> Float64(x, RoundDown), values),
    map(x -> Float64(x, RoundUp), values),
)
function ads_product_bounds(lo, hi, xlo, xhi)
    lo == hi == 0 && return (0.0, 0.0)
    lo == hi == 1 && return (xlo, xhi)
    lo == hi == -1 && return (-xhi, -xlo)
    products = (lo * xlo, lo * xhi, hi * xlo, hi * xhi)
    any(isnan, products) && return (-Inf, Inf)
    return prevfloat(minimum(products)), nextfloat(maximum(products))
end
function ads_sum_bounds(lo, hi, xlo, xhi)
    lo == hi == 0 && return (xlo, xhi)
    xlo == xhi == 0 && return (lo, hi)
    lower, upper = lo + xlo, hi + xhi
    (isnan(lower) || isnan(upper)) && return (-Inf, Inf)
    return prevfloat(lower), nextfloat(upper)
end
function ads_project_bounds(matrix_bounds, point_bounds)
    lower, upper = matrix_bounds
    xlo, xhi = point_bounds
    lo, hi = zeros(size(lower, 1)), zeros(size(lower, 1))
    for i in axes(lower, 1), j in axes(lower, 2)
        a, b = ads_product_bounds(lower[i, j], upper[i, j], xlo[j], xhi[j])
        lo[i], hi[i] = ads_sum_bounds(lo[i], hi[i], a, b)
    end
    return lo, hi
end

polygon_cross(a, b, c) = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
polygon_signed_area(v) = sum(v[i][1] * v[mod1(i + 1, length(v))][2] - v[i][2] * v[mod1(i + 1, length(v))][1] for i in eachindex(v)) / 2

"""
    ConvexPolygon(vertices)

A nondegenerate convex polygon with ordered 2D vertices, clockwise or
counterclockwise. Geometry uses exact rationals; floating inputs denote their
stored binary values. Rejects nonconvex, self-intersecting and invalid domains.
Construction owns its storage. Access vertices with `polygon_vertices`.
"""
struct ConvexPolygon{T <: AbstractFloat}
    _vertices::Tuple{Vararg{PolygonPoint}}
    ConvexPolygon{T}(v, ::Val{:owned}) where {T <: AbstractFloat} = new{T}(v)
end
function ConvexPolygon(vertices)
    vertices isa Union{Tuple, AbstractVector} && length(vertices) >= 3 || throw(ArgumentError("Supply at least three ordered polygon vertices"))
    all(v -> v isa Union{Tuple, AbstractVector} && length(v) == 2, vertices) || throw(DimensionMismatch("Polygon vertices require two coordinates"))
    v = Tuple((polygon_real(p[1]), polygon_real(p[2])) for p in vertices)
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten(vertices))
    T <: AbstractFloat && isconcretetype(T) || throw(ArgumentError("Unsupported polygon coordinate type"))
    length(unique(v)) == length(v) || throw(ArgumentError("Repeated polygon vertices"))
    area = polygon_signed_area(v)
    iszero(area) && throw(ArgumentError("Polygon must have positive area"))
    area < 0 && (v = reverse(v))
    all(i -> all(p -> polygon_cross(v[i], v[mod1(i + 1, length(v))], p) >= 0, v), eachindex(v)) ||
        throw(ArgumentError("Require ordered vertices of a convex, simple polygon"))
    return ConvexPolygon{T}(deepcopy(v), Val(:owned))
end
"""
    polygon_vertices(polygon::ConvexPolygon)

Return independent exact rational vertices in counterclockwise order.
"""
polygon_vertices(p::ConvexPolygon) = collect(deepcopy(p._vertices))
"""
    domain_area(polygon::ConvexPolygon)

Return the exact rational area of the physical polygon.
"""
domain_area(p::ConvexPolygon) = polygon_signed_area(p._vertices)
Base.copy(p::ConvexPolygon{T}) where {T} = ConvexPolygon{T}(deepcopy(p._vertices), Val(:owned))
Base.deepcopy_internal(p::ConvexPolygon, copies::IdDict) = get!(() -> copy(p), copies, p)
Base.show(io::IO, p::ConvexPolygon) = print(io, "ConvexPolygon(", length(p._vertices), " vertices)")
polygon_contains(p::ConvexPolygon, x) = all(i -> polygon_cross(p._vertices[i], p._vertices[mod1(i + 1, length(p._vertices))], x) >= 0, eachindex(p._vertices))

function polygon_filter(p::ConvexPolygon)
    planes = Matrix{PolygonReal}(undef, length(p._vertices), 3)
    for i in eachindex(p._vertices)
        a, b = p._vertices[i], p._vertices[mod1(i + 1, length(p._vertices))]
        dx, dy = b[1] - a[1], b[2] - a[2]
        planes[i, :] .= (dy, -dx, dy * a[1] - dx * a[2])
    end
    return ads_float_bounds(planes)
end
function polygon_contains_filtered(p::ConvexPolygon, point, bounds, planes)
    lower, upper = planes
    lo, hi = ads_project_bounds((view(lower, :, 1:2), view(upper, :, 1:2)), bounds)
    for i in eachindex(lo)
        hi[i] <= lower[i, 3] && continue
        lo[i] > upper[i, 3] && return false
        a, b = p._vertices[i], p._vertices[mod1(i + 1, length(p._vertices))]
        polygon_cross(a, b, polygon_real.(point)) >= 0 || return false
    end
    return true
end

function polygon_box(lower, upper)
    length(lower) == length(upper) == 2 || throw(DimensionMismatch("Polygonal ADS requires two coordinate bounds"))
    a, b = polygon_real.(lower), polygon_real.(upper)
    all(a .< b) || throw(ArgumentError("Polygonal ADS requires positive area; use box ADS for fixed coordinates"))
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    v = ((a[1], a[2]), (b[1], a[2]), (b[1], b[2]), (a[1], b[2]))
    return ConvexPolygon{T}(deepcopy(v), Val(:owned))
end

# Clip against a*x+b*y <= cut. Every intersection is exact; shared faces
# belong to both closed children and their union equals the closed parent.
function polygon_clip_points(v, normal, cut)
    isempty(v) && return ()
    result = PolygonPoint[]
    previous = last(v)
    previous_value = normal[1] * previous[1] + normal[2] * previous[2] - cut
    for current in v
        value = normal[1] * current[1] + normal[2] * current[2] - cut
        if (value <= 0) != (previous_value <= 0)
            t = previous_value / (previous_value - value)
            crossing = (previous[1] + t * (current[1] - previous[1]), previous[2] + t * (current[2] - previous[2]))
            (isempty(result) || last(result) != crossing) && push!(result, crossing)
        end
        value <= 0 && (isempty(result) || last(result) != current) && push!(result, current)
        previous, previous_value = current, value
    end
    length(result) > 1 && first(result) == last(result) && pop!(result)
    return Tuple(result)
end
function polygon_clip(p::ConvexPolygon{T}, normal, cut) where {T}
    result = polygon_clip_points(p._vertices, normal, cut)
    length(result) >= 3 && polygon_signed_area(result) > 0 || return nothing
    return ConvexPolygon{T}(result, Val(:owned))
end
function polygon_intersection(a::ConvexPolygon, b::ConvexPolygon)
    result = a
    v = b._vertices
    for i in eachindex(v)
        x, y = v[i], v[mod1(i + 1, length(v))]
        normal = (y[2] - x[2], x[1] - y[1])
        result = polygon_clip(result, normal, normal[1] * x[1] + normal[2] * x[2])
        result === nothing && return nothing
    end
    return result
end

function polygon_frame(directions)
    directions isa AbstractMatrix && size(directions) == (2, 2) || throw(DimensionMismatch("directions must be a 2x2 matrix of projection rows"))
    B = polygon_real.(directions)
    determinant = B[1, 1] * B[2, 2] - B[1, 2] * B[2, 1]
    iszero(determinant) && throw(ArgumentError("Splitting directions must be linearly independent"))
    A = [B[2, 2] -B[1, 2]; -B[2, 1] B[1, 1]] ./ determinant
    return B, A
end
polygon_project(B, x) = (B[1, 1] * x[1] + B[1, 2] * x[2], B[2, 1] * x[1] + B[2, 2] * x[2])
function polygon_projection_bounds(p, B)
    points = map(x -> polygon_project(B, x), p._vertices)
    return ntuple(i -> minimum(x -> x[i], points), 2), ntuple(i -> maximum(x -> x[i], points), 2)
end
