module PolygonQueryTestSupport
using DifferentialAlgebra
const DA = DifferentialAlgebra
const Exact = Rational{BigInt}

# The previous exact traversal is an oracle for left-biased point routing.
function exact_patch(a, point)
    x = Exact.(DA.continuity_primal.(point))
    node = first(a.nodes)
    while node.patch == 0
        value = sum(node.normal[i] * x[i] for i in 1:2)
        node = a.nodes[value <= node.midpoint ? node.left : node.right]
    end
    return node.patch
end
function exact_evaluate(a, point)
    p = DA.ads_patches(a)[exact_patch(a, point)]._patch
    T = promote_type(eltype(p.center), eltype(point))
    projected = T.(a._directions) * T[point...]
    normalized = [
        iszero(p.radius[i]) ? zero(T) : (projected[i] - p.center[i]) / p.radius[i]
            for i in 1:2
    ]
    result = evaluate(p.map, normalized)
    return a._scalar ? only(result) : result
end
function face_points(a)
    points = DA.PolygonPoint[]
    for p in DA.ads_patches(a)
        v = polygon_vertices(domain(p))
        append!(points, v)
        for i in eachindex(v)
            x, y = v[i], v[mod1(i + 1, length(v))]
            push!(points, ((x[1] + y[1]) / 2, (x[2] + y[2]) / 2))
        end
    end
    return unique(points)
end
end
