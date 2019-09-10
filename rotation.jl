"""
    rotation_matrix_euler!(mat, phi, theta, psi, inv=false)

Overwrite `mat[1:3, 1:3]` with the rotation matrix specified by the
euler angles `phi`, `theta` and `psi`.
"""
function rotation_matrix_euler!(mat, phi, theta, psi, inv=false)
    # convert degrees to radians
    alpha, beta, gamma = map(x -> x * pi / 180, [phi, theta, psi])

    c1, c2, c3 = cos(alpha), cos(beta), cos(gamma)
    s1, s2, s3 = sin(alpha), sin(beta), sin(gamma)

    # define rotation matrix (zyz convention)
    if inv
        mat[1, 1] = c1 * c2 * c3 - s1 * s3
        mat[1, 2] = c1 * s3 + c2 * c3 * s1
        mat[1, 3] = -c3 * s2
        mat[2, 1] = -c3 * s1 - c1 * c2 * s3
        mat[2, 2] = c1 * c3 - c2 * s1 * s3
        mat[2, 3] = s2 * s3
        mat[3, 1] = c1 * s2
        mat[3, 2] = s1 * s2
        mat[3, 3] = c2
    else
        mat[1, 1] = c1 * c2 * c3 - s1 * s3
        mat[1, 2] = -c3 * s1 - c1 * c2 * s3
        mat[1, 3] = c1 * s2
        mat[2, 1] = c1 * s3 + c2 * c3 * s1
        mat[2, 2] = c1 * c3 - c2 * s1 * s3
        mat[2, 3] = s1 * s2
        mat[3, 1] = -c3 * s2
        mat[3, 2] = s2 * s3
        mat[3, 3] = c2
    end
end
