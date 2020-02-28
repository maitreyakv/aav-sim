/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_QUADCOPTERDYNAMICS_H
#define MPCSIM_QUADCOPTERDYNAMICS_H


#include "Dynamics.h"


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

    // Product of k and l
    double m_kl;

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
    {
        // Compute product of k and l
        this->m_kl = k * l;
    }

    /**
     * Computes the state derivative dxdt = F(x,u) of the Dynamics
     *
     * @param x Vector with the state "x" to evaluate F(x,u) at
     * @param u Vector with the control input "u" to evaluate F(x,u) at
     * @param t The time to evaluate F(x,u) at (not used since F is time-invariant)
     * @return Vector containing the state derivative dx/dt = F(x,u)
     */
    arma::vec F(arma::vec x, arma::vec u, double t);

    /**
     * Computes the gradient of the state derivative F_x of the Dynamics with respect to the state,
     * returned in the form
     *
     *      Phi = I + F_x(x,u) * dt
     *
     * @param x Vector with the state "x" to evaluate F_x(x,u) at
     * @param u Vector with the control input "u" to evaluate F_x(x,u) at
     * @param dt The discretized time step in the DDP algorithm, t(k+1) - t(k)
     * @return Matrix containing Phi
     */
    arma::mat Phi(arma::vec x, arma::vec u, double dt);

    /**
     * Computes the gradient of the state derivative F_x of the Dynamics with respect to the control,
     * returned in the form
     *
     *      Beta = F_u(x,u) * dt
     *
     * @param x Vector with the state "x" to evaluate F_x(x,u) at
     * @param u Vector with the control input "u" to evaluate F_x(x,u) at
     * @param dt The discretized time step in the DDP algorithm, t(k+1) - t(k)
     * @return Matrix containing Beta
     */
    arma::mat Beta(arma::vec x, arma::vec u, double dt);

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


#endif //MPCSIM_QUADCOPTERDYNAMICS_H
