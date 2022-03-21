from casadi import SX, cross, dot, norm_2, vertcat


def sx_quat_normalize(q):

    return q / norm_2(q)


def sx_quat_multiply(q, p):

    s1 = q[0]
    v1 = q[1:4]
    s2 = p[0]
    v2 = p[1:4]
    s = s1 * s2 - dot(v1, v2)
    v = s1 * v2 + s2 * v1 + cross(v1, v2)
    return vertcat(s, v)


def sx_quat_inverse(q):

    return SX([1, -1, -1, -1]) * q / (norm_2(q)**2)


def quat_derivative(q, w):

    q = sx_quat_normalize(q)
    return sx_quat_multiply(q, vertcat(SX(0), w)) / 2


def rotate_vec3(v, q):

    q = sx_quat_normalize(q)
    p = vertcat(SX(0), v)

    p_rotated = sx_quat_multiply(sx_quat_multiply(q, p), sx_quat_inverse(q))
    return p_rotated[1:4]


def quat_err(q, p):

    # assume q = q_delta * p, q_delta = [w, q_err] = q * p^(-1)
    q = sx_quat_normalize(q)
    p = sx_quat_normalize(p)
    q_d = sx_quat_multiply(q, sx_quat_inverse(p))
    return q_d[1:4]
