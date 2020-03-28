/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAVSIM_QUADCOPTERDYNAMICS_H
#define AAVSIM_QUADCOPTERDYNAMICS_H


#include "../simulation/Dynamics.h"

// TODO: Update documentation

class QuadcopterDynamics : public Dynamics {

private:
    // Mass of quadcopter
    double m_m;

    // Gravitation acceleration (typically 9.8 m/s^2)
    double m_g;

    // Distance from a propeller to center of quadcopter
    double m_l;

    // Gain of the motor-propeller, multiplied square of angular speed to get thrust
    double m_k;

    // Friction parameter of the motors
    double m_b;

    // Components of principle moment of inertia
    double m_I_xx;
    double m_I_yy;
    double m_I_zz;

public:
    QuadcopterDynamics(double m, double g, double l, double k, double b, double I_xx, double I_yy, double I_zz)
        : m_m(m),
          m_g(g),
          m_l(l),
          m_k(k),
          m_b(b),
          m_I_xx(I_xx),
          m_I_yy(I_yy),
          m_I_zz(I_zz)
    {}

    /**
     * Computes the state derivative dxdt = F(x,u) of the Dynamics
     *
     * @param x Vector with the state "x" to evaluate F(x,u) at
     * @param u Vector with the control input "u" to evaluate F(x,u) at
     * @param t The time to evaluate F(x,u) at (not used since F is time-invariant)
     * @return Vector containing the state derivative dx/dt = F(x,u)
     */
    arma::vec F(arma::vec x, arma::vec u, double t);

    arma::vec f(arma::vec x, arma::vec u, double t, double dt);

    arma::mat f_x(arma::vec x, arma::vec u, double t, double dt);

    arma::mat f_u(arma::vec x, arma::vec u, double t, double dt);

    arma::cube f_xx(arma::vec x, arma::vec u, double t, double dt);

    arma::cube f_uu(arma::vec x, arma::vec u, double t, double dt);

    arma::cube f_ux(arma::vec x, arma::vec u, double t, double dt);

    /**
     * Returns the number of elements (dimensionality) of the state vector, which is 12 for this system:
     *
     *  x = [X, Y, Z, X_dot, Y_dot, Z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]^T
     *
     * @return 4
     */
    int getStateDimension() {return 12;}

    /**
     * Returns the number of elements (dimensionality) of the control vector, which is 4 for this system:
     *
     *  u = [omega_1, omega_2, omega_3, omega_4]^T
     *
     * @return 4
     */
    int getControlDimension() {return 4;}
};


#endif //AAVSIM_QUADCOPTERDYNAMICS_H
