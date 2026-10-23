# Planar circular restricted three-body model shared by the propagation and
# uncertainty tutorials. Distance, total mass and primary angular speed are one;
# the state is (x, y, vx, vy) in the rotating frame.
function cr3bp(u, μ, t)
    x, y, vx, vy = u
    dx1, dx2 = x + μ, x - 1 + μ
    r13 = (dx1^2 + y^2)^(3 // 2)
    r23 = (dx2^2 + y^2)^(3 // 2)
    return [
        vx, vy,
        2vy + x - (1 - μ) * dx1 / r13 - μ * dx2 / r23,
        -2vx + y - (1 - μ) * y / r13 - μ * y / r23,
    ]
end

# The conserved Jacobi constant is -2 times the rotating-frame energy.
function jacobi(u, μ)
    x, y, vx, vy = u
    r1 = sqrt((x + μ)^2 + y^2)
    r2 = sqrt((x - 1 + μ)^2 + y^2)
    return x^2 + y^2 + 2 * ((1 - μ) / r1 + μ / r2) - vx^2 - vy^2
end
