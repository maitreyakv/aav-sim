from sympy import *
from sympy.simplify.cse_main import cse
import shutil
import numpy as np

def insert_code(tag, code, last):
    with open("QuadcopterDynamics.cpp") as f:
        string = ""
        # Subexpression
        if tag == "@SUBEXPR":
            if last:
                string = "double %s = %s;" % (str(code[0]), ccode(code[1]))
            else:
                string = "double %s = %s;\n    @SUBEXPR" % (str(code[0]), ccode(code[1]))
            new_source = f.read().replace(tag, string, 1)
        # Tensor component
        elif "f_xx" in tag or "f_uu" in tag or "f_ux" in tag:
            if last:
                string = "%s = %s;" % (tag, ccode(code))
            else:
                string = "%s = %s;\n    @TENSOR" % (tag, ccode(code))
            new_source = f.read().replace("@TENSOR", string, 1)
        # Jacobian component
        elif "f_x" in tag or "f_u" in tag:
            if last:
                string = "%s = %s;" % (tag, ccode(code))
            else:
                string = "%s = %s;\n    @JACOBIAN" % (tag, ccode(code))
            new_source = f.read().replace("@JACOBIAN", string, 1)
        # Function f component
        else:
            string = str(ccode(code) + ";")
            new_source = f.read().replace(tag, string, 1)

    with open("QuadcopterDynamics.cpp", "w") as f:
        f.write(new_source)

if __name__ == "__main__":
    # Create CPP files
    shutil.copy("QuadcopterDynamics.cpp.in", "QuadcopterDynamics.cpp")

    # Define independent variables (state variable and parameters)
    X, Y, Z, phi, theta, psi = symbols('x[0] x[1] x[2] x[6] x[7] x[8]')
    X_dot, Y_dot, Z_dot, phi_dot, theta_dot, psi_dot = \
        symbols('x[3] x[4] x[5] x[9] x[10] x[11]')

    x = Matrix([X, Y, Z, X_dot, Y_dot, Z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot])

    # Define quadcopter parameters
    m, g, l, k, b, I_xx, I_yy, I_zz = symbols('this->m_m this->m_g this->m_l this->m_k this->m_b \
                                               this->m_I_xx this->m_I_yy this->m_I_zz')

    # Define control inputs
    u_1, u_2, u_3, u_4 = symbols('u[0] u[1] u[2] u[3]')

    u = Matrix([u_1, u_2, u_3, u_4])

    # Define discretization time step
    dt = symbols('dt')

    gamma_1 = (sqrt(m*g/(4*k)) + u_1)**2
    gamma_2 = (sqrt(m*g/(4*k)) + u_2)**2
    gamma_3 = (sqrt(m*g/(4*k)) + u_3)**2
    gamma_4 = (sqrt(m*g/(4*k)) + u_4)**2

    # Intermediate expressions
    s_phi = sin(phi)
    s_theta = sin(theta)
    s_psi = sin(psi)
    c_phi = cos(phi)
    c_theta = cos(theta)
    c_psi = cos(psi)

    A = Matrix([[1, 0, -s_theta],
                [0, c_phi, c_theta*s_phi],
                [0, -s_phi, c_theta*c_phi]])

    omega = A * Matrix([phi_dot, theta_dot, psi_dot])

    R = Matrix([[c_phi*c_theta, c_phi*s_theta*s_psi - c_psi*s_phi, s_phi*s_psi + c_phi*c_psi*s_theta],
                [c_theta*s_phi, c_phi*c_psi + s_phi*s_theta*s_psi, c_psi*s_phi*s_theta - c_phi*s_psi],
                [-s_theta, c_theta*s_psi, c_theta*c_psi]])

    T = Matrix([0, 0, k*(gamma_1 + gamma_2 + gamma_3 + gamma_4)])

    tau = Matrix([l*k*(gamma_1 - gamma_3),
                  l*k*(gamma_2 - gamma_4),
                  b*(gamma_1 - gamma_2 + gamma_3 - gamma_4)])

    x_ddot = (Matrix([0, 0, -m*g]) + R*T) / m

    I = Matrix([[I_xx, 0, 0], [0, I_yy, 0], [0, 0, I_zz]])

    omega_dot = I.inv() * (tau - omega.cross(I*omega))

    Theta_dot = Matrix([phi_dot, theta_dot, psi_dot])

    Theta_ddot = A.inv()*(omega_dot - \
        Matrix([[0, 0, -c_theta*theta_dot],
                [0, -s_phi*phi_dot, -s_theta*theta_dot*s_phi + c_theta*c_phi*phi_dot],
                [0, -c_phi*phi_dot, -s_theta*theta_dot*c_phi - c_theta*s_phi*phi_dot]])
                    *Theta_dot)

    F = Matrix([X_dot, Y_dot, Z_dot, \
                x_ddot[0], x_ddot[1], x_ddot[2], \
                phi_dot, theta_dot, psi_dot, \
                Theta_ddot[0], Theta_ddot[1], Theta_ddot[2]])

    f = x + F * dt

    f_x = f.jacobian(x)
    f_u = f.jacobian(u)

    f_xx = np.zeros((12, 12, 12), dtype=Symbol)
    for i in range(12):
        for j in range(12):
            for k in range(12):
                f_xx[i,j,k] = diff(f_x[i,j], x[k])

    f_uu = np.zeros((12, 4, 4), dtype=Symbol)
    for i in range(12):
        for j in range(4):
            for k in range(4):
                f_uu[i,j,k] = diff(f_u[i,j], u[k])

    f_ux = np.zeros((12, 4, 12), dtype=Symbol)
    for i in range(12):
        for j in range(4):
            for k in range(12):
                f_ux[i,j,k] = diff(f_u[i,j], x[k])

    # Write function F
    subexp, F_reduced = cse(F, optimizations='basic')
    F_reduced = F_reduced[0]

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        insert_code("@F_%d" % i, F_reduced[i], True)

    # Write function f
    subexp, f_reduced = cse(f, optimizations='basic')
    f_reduced = f_reduced[0]

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        insert_code("@f_%d" % i, f_reduced[i], True)

    # Write function f_x
    subexp, f_x_reduced = cse(f_x, optimizations='basic')
    f_x_reduced = f_x_reduced[0]

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        for j in range(12):
            tag = "f_x(%d,%d)" % (i,j)
            if not f_x_reduced[i,j] == 0:
                insert_code(tag, f_x_reduced[i,j], i == 11 and j == 11)

    # Write function f_u
    subexp, f_u_reduced = cse(f_u, optimizations='basic')
    f_u_reduced = f_u_reduced[0]

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        for j in range(4):
            tag = "f_u(%d,%d)" % (i,j)
            if not f_u_reduced[i,j] == 0:
                insert_code(tag, f_u_reduced[i,j], i == 11 and j == 3)

    # Write function f_xx
    subexp, f_xx_reduced = cse(f_xx.flatten().tolist(), optimizations='basic')

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        for j in range(12):
            for k in range(12):
                tag = "f_xx(%d,%d,%d)" % (i,j,k)
                if not f_xx_reduced[144*i + 12*j + k] == 0:
                    insert_code(tag, f_xx_reduced[144*i + 12*j + k], i == 11 and j == 11 and k == 11)

    # Write function f_uu
    subexp, f_uu_reduced = cse(f_uu.flatten().tolist(), optimizations='basic')

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        for j in range(4):
            for k in range(4):
                tag = "f_uu(%d,%d,%d)" % (i,j,k)
                if not f_uu_reduced[16*i + 4*j + k] == 0:
                    insert_code(tag, f_uu_reduced[16*i + 4*j + k], i == 11 and j == 3 and k == 3)

    # Write function f_ux
    subexp, f_ux_reduced = cse(f_ux.flatten().tolist(), optimizations='basic')

    for exp in subexp:
        insert_code("@SUBEXPR", exp, exp == subexp[-1])

    for i in range(12):
        for j in range(4):
            for k in range(12):
                tag = "f_ux(%d,%d,%d)" % (i,j,k)
                if not f_ux_reduced[48*i + 12*j + k] == 0:
                    insert_code(tag, f_ux_reduced[48*i + 12*j + k], i == 11 and j == 3 and k == 11)

    # Remove all remaining tags
    with open("QuadcopterDynamics.cpp") as f:
        lines = f.readlines()

    with open("QuadcopterDynamics.cpp", "w") as f:
        for line in lines:
            if not "@SUBEXPR" in line and not "@JACOBIAN" in line and not "@TENSOR" in line:
                f.write(line)
