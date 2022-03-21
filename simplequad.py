import math
import shutil
import timeit
from math import cos, pi, sin
from os import environ, path, removedirs, system

import matplotlib.pyplot as plt
import numpy as np
from acados_template import (AcadosModel, AcadosOcp, AcadosOcpSolver,
                             AcadosSim, AcadosSimSolver)
from casadi import SX, Function, cross, diag, dot, inv, norm_2, sqrt, vertcat

from sx_quaternion import (quat_derivative, quat_err, rotate_vec3,
                           sx_quat_inverse, sx_quat_multiply,
                           sx_quat_normalize)


def get_simplequad_model():

    m = 1
    g = SX([0, 0, 10])
    J_xx = 1
    J_xy = 0
    J_xz = 0
    J_yy = 1
    J_yz = 0
    J_zz = 1
    J = SX([
        [J_xx, J_xy, J_xz],
        [J_xy, J_yy, J_yz],
        [J_xz, J_yz, J_zz]
    ])
    J_inv = inv(J)

    p = SX.sym('p', 3)  # position in world NED frame (m)
    q = SX.sym('q', 4)  # attitude expressed in quaternion
    v = SX.sym('v', 3)  # velocity in world frame (m/s)
    w = SX.sym('w', 3)  # angular velocity in body FRD frame (rad/s)
    f = SX.sym('f', 3)  # actuator force in body frame (N)
    t = SX.sym('t', 3)  # actuator torque in body frame (Nm)

    p_ref = SX.sym('p_ref', 3)
    q_ref = SX.sym('q_ref', 4)
    v_ref = SX.sym('v_ref', 3)
    w_ref = SX.sym('w_ref', 3)

    p_dot = SX.sym('p_dot', 3)
    q_dot = SX.sym('q_dot', 4)
    v_dot = SX.sym('v_dot', 3)
    w_dot = SX.sym('w_dot', 3)

    x = vertcat(p, q, v, w)
    x_dot = vertcat(p_dot, q_dot, v_dot, w_dot)
    u = vertcat(f, t)
    p = vertcat(p_ref, q_ref, v_ref, w_ref)

    f_expl = vertcat(
        v,
        quat_derivative(q, w),
        rotate_vec3(f, q) / m + g,
        J_inv @ (t - cross(w, J @ w)),
    )
    f_impl = x_dot - f_expl

    model = AcadosModel()
    model.f_expl_expr = f_expl
    model.f_impl_expr = f_impl
    model.x = x
    model.xdot = x_dot
    model.u = u
    model.p = p
    model.name = 'simplequad'

    print("model summary")
    print('m =', m)
    print('J =', J)
    print('x =', x)
    print('u =', u)
    print('p =', p)

    return model

def e(x, ref):

    p_ref = ref[0:3]
    q_ref = ref[3:7]
    v_ref = ref[7:10]
    w_ref = ref[10:13]

    p_err = x[0:3] - p_ref
    q_err = quat_err(x[3:7], q_ref)
    v_err = x[7:10] - v_ref
    w_err = x[10:13] - w_ref

    return vertcat(p_err, q_err, v_err, w_err)

def formulate_simplequad_ocp():

    Q_p = SX([1, 1, 1])
    Q_q = SX([1, 1, 1])
    Q_v = SX([1, 1, 1])
    Q_w = SX([1, 1, 1])
    Q = diag(vertcat(Q_p, Q_q, Q_v, Q_w))

    Q_e = Q

    R_f = SX([1, 1, 1])
    R_t = SX([1, 1, 1])
    R = diag(vertcat(R_f, R_t))

    u_lb = np.array([0, 0, -40, -20, -20, -20])
    u_ub = np.array([0, 0, 0, 20, 20, 20])
    x_lb = np.array([-10, -10, -10, -60, -60, -60])
    x_ub = np.array([10, 10, 10, 60, 60, 60])

    ocp = AcadosOcp()
    ocp.model = get_simplequad_model()
    
    ocp.dims.N = 20
    x = ocp.model.x
    u = ocp.model.u
    p = ocp.model.p

    ocp.model.cost_expr_ext_cost = e(x, p).T @ Q @ e(x, p) + u.T @ R @ u
    ocp.model.cost_expr_ext_cost_e = e(x, p).T @ Q_e @ e(x, p)

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    ocp.constraints.constr_type = 'BGH'
    ocp.constraints.idxbu = np.array([0, 1, 2, 3, 4, 5])
    ocp.constraints.lbu = u_lb
    ocp.constraints.ubu = u_ub
    ocp.constraints.idxbx = np.array([7, 8, 9, 10, 11, 12])
    ocp.constraints.lbx = x_lb
    ocp.constraints.ubx = x_ub
    # ocp.constraints.x0 = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    ocp.solver_options.tf = 1.0
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp.solver_options.qp_solver_iter_max = 100
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_HPIPM'
    # ocp.solver_options.nlp_solver_type = 'SQP'
    # ocp.solver_options.levenberg_marquardt = 1e-3
    # ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'

    ocp.acados_include_path = environ['ACADOS_SOURCE_DIR'] + '/include'
    ocp.acados_lib_path = environ['ACADOS_SOURCE_DIR'] + '/lib'
    ocp.code_export_directory = './acados_export'

    ocp.parameter_values = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    return ocp


def get_sim_integrator(ocp, dt):

    sim = AcadosSim()

    sim.model = ocp.model
    sim.acados_include_path = ocp.acados_include_path
    sim.acados_lib_path = ocp.acados_lib_path
    sim.parameter_values = ocp.parameter_values
    sim.code_export_directory = ocp.code_export_directory

    sim.solver_options.T = dt
    sim.solver_options.newton_iter = 3
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 1

    return sim


if __name__ == '__main__':

    np.set_printoptions(linewidth=1e6, suppress=True)

    # ocp
    ocp = formulate_simplequad_ocp()

    if path.exists(ocp.code_export_directory):
        shutil.rmtree(ocp.code_export_directory)

    print('================= generate ocp solver c code =================')
    ocp_json_path = ocp.code_export_directory + '/acados_ocp_' + ocp.model.name + '.json'
    ocpSolver = AcadosOcpSolver(ocp, ocp_json_path, build=True)

    # sim
    control_rate = 100
    dt = 1 / control_rate
    sim = get_sim_integrator(ocp, dt)

    print('================= generate sim solver c code =================')
    sim_json_path = sim.code_export_directory + '/acados_sim_' + sim.model.name + '.json'
    simSolver = AcadosSimSolver(sim, sim_json_path, build=True)

    # test sim
    print('================= test sim solver =================')
    x = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    u = np.array([0, 0, -11, 0, 0, 1])
    
    for i in range(1001):
        t = i * 0.01
        if i == 1000 or i == 100 or i == 10:
            print("t =", t)
            print(x)
        simSolver.set('x', x)
        simSolver.set('u', u)
        simSolver.solve()
        x = simSolver.get('x')

    # test ocp
    print('================= test ocp solver =================')
    x = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    print('x =', x)

    y = np.zeros([13, 21])
    y[3, :] = 1

    for i in range(21):
        t = i * 0.05
        y[2, i] = -0.5 + 0.5 * cos(pi * t)
        y[3, i] = 1
        y[9, i] = -0.5 * pi * sin(pi * t)
        ocpSolver.set(i, 'p', y[:, i])

    print('ref =')
    print(y)

    ocpSolver.set(0, 'x', x)
    
    ocp_status = ocpSolver.solve()
    ocpSolver.print_statistics()

    if ocp_status != 0:
        print('OCP_SOL: ocp_status =', ocp_status, ', exiting')
        quit()

    u = ocpSolver.get(0, 'u')
    print(u)


    
