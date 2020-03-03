from sympy import *
from sympy.simplify.cse_main import cse
import shutil

def convert_to_cpp_code(string):
    string = string.replace("X_dot", "x[3]")
    string = string.replace("Y_dot", "x[4]")
    string = string.replace("Z_dot", "x[5]")

    string = string.replace("X", "x[0]")
    string = string.replace("Y", "x[1]")
    string = string.replace("Z", "x[2]")

    string = string.replace("phi_dot", "x[9]")
    string = string.replace("theta_dot", "x[10]")
    string = string.replace("psi_dot", "x[11]")

    string = string.replace("phi", "x[6]")
    string = string.replace("theta", "x[7]")
    string = string.replace("psi", "x[8]")

    string = string.replace("u_1", "u[0]")
    string = string.replace("u_2", "u[1]")
    string = string.replace("u_3", "u[2]")
    string = string.replace("u_4", "u[3]")

    string = string.replace("phi_dot", "x[9]")
    string = string.replace("theta_dot", "x[10]")
    string = string.replace("psi_dot", "x[11]")

    string = string.replace("m", "this->m_m")
    string = string.replace("b", "this->m_b")
    string = string.replace("g", "this->m_g")
    string = string.replace("I_xx", "this->m_I_xx")
    string = string.replace("I_yy", "this->m_I_yy")
    string = string.replace("I_zz", "this->m_I_zz")

    string = string.replace("k", "this->m_k")
    string = string.replace("l", "this->m_l")

    for n in range(100):
        string = string.replace("x%d**2" % n, "pow(x%d, 2)" % n)

    return string

def insert_code(tag, code, last):
    with open("QuadcopterDynamics.cpp") as f:
        string = ""
        if tag == "@SUBEXPR":
            if last:
                string = "double %s = %s;" % (str(code[0]), convert_to_cpp_code(str(code[1])))
            else:
                string = "double %s = %s;\n    @SUBEXPR" % (str(code[0]), convert_to_cpp_code(str(code[1])))
            new_source = f.read().replace(tag, string, 1)
        elif "F_" in tag:
            if last:
                string = "%s = %s;" % (tag[1:], convert_to_cpp_code(str(code)))
            else:
                string = "%s = %s;\n    @JACOBIAN" % (tag[1:], convert_to_cpp_code(str(code)))
            new_source = f.read().replace("@JACOBIAN", string, 1)
        else:
            string = str(convert_to_cpp_code(str(code)) + ";")
            new_source = f.read().replace(tag, string, 1)

    with open("QuadcopterDynamics.cpp", "w") as f:
        f.write(new_source)

if __name__ == "__main__":
    # Create CPP files
    shutil.copy("QuadcopterDynamics.cpp.in", "QuadcopterDynamics.cpp")

    # Define independent variables (state variable and parameters)
    X, Y, Z, phi, theta, psi = symbols('X Y Z phi theta psi')
    X_dot, Y_dot, Z_dot, phi_dot, theta_dot, psi_dot = \
        symbols('X_dot Y_dot Z_dot phi_dot theta_dot psi_dot')

    m, g, l, k, b, I_xx, I_yy, I_zz = symbols('m g l k b I_xx I_yy I_zz')
    u_1, u_2, u_3, u_4 = symbols('u_1 u_2 u_3 u_4')

    gamma_1 = 1000*u_1 + m*g/(4*k)
    gamma_2 = 1000*u_2 + m*g/(4*k)
    gamma_3 = 1000*u_3 + m*g/(4*k)
    gamma_4 = 1000*u_4 + m*g/(4*k)

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

    F_x = F.jacobian(Matrix([X, Y, Z, X_dot, Y_dot, Z_dot, \
                             phi, theta, psi, phi_dot, theta_dot, psi_dot]))

    F_u = F.jacobian(Matrix([u_1, u_2, u_3, u_4]))

    # Write function F
    eqns = {}

    eqns["@X_DDOT"] = x_ddot[0]
    eqns["@Y_DDOT"] = x_ddot[1]
    eqns["@Z_DDOT"] = x_ddot[2]

    eqns["@PHI_DDOT"] = Theta_ddot[0]
    eqns["@THETA_DDOT"] = Theta_ddot[1]
    eqns["@PSI_DDOT"] = Theta_ddot[2]

    subexp, eqns_reduced = cse(eqns.values(), optimizations='basic')
    i = 0
    for tag in eqns.keys():
        eqns[tag] = eqns_reduced[i]
        i += 1

    for exp in subexp:
        if exp == subexp[-1]:
            insert_code("@SUBEXPR", exp, True)
        else:
            insert_code("@SUBEXPR", exp, False)

    for tag, eqn in eqns.items():
        insert_code(tag, eqn, True)

    # Write function Phi
    eqns = {}

    for i in range(12):
        for j in range(12):
            eqns["@F_x(%d,%d)" % (i,j)] = F_x[i,j]

    subexp, eqns_reduced = cse(eqns.values(), optimizations='basic')
    i = 0
    for tag in eqns.keys():
        eqns[tag] = eqns_reduced[i]
        i += 1

    for exp in subexp:
        if exp == subexp[-1]:
            insert_code("@SUBEXPR", exp, True)
        else:
            insert_code("@SUBEXPR", exp, False)

    for i in range(12):
        for j in range(12):
            tag = "@F_x(%d,%d)" % (i,j)
            if not eqns[tag] == 0:
                if i == 11 and j == 11:
                    insert_code(tag, eqns[tag], True)
                else:
                    insert_code(tag, eqns[tag], False)

    # Write function Beta
    eqns = {}

    for i in range(12):
        for j in range(4):
            eqns["@F_u(%d,%d)" % (i,j)] = F_u[i,j]

    subexp, eqns_reduced = cse(eqns.values(), optimizations='basic')
    i = 0
    for tag in eqns.keys():
        eqns[tag] = eqns_reduced[i]
        i += 1

    for exp in subexp:
        if exp == subexp[-1]:
            insert_code("@SUBEXPR", exp, True)
        else:
            insert_code("@SUBEXPR", exp, False)

    for i in range(12):
        for j in range(4):
            tag = "@F_u(%d,%d)" % (i,j)
            if not eqns[tag] == 0:
                if i == 11 and j == 3:
                    insert_code(tag, eqns[tag], True)
                else:
                    insert_code(tag, eqns[tag], False)
