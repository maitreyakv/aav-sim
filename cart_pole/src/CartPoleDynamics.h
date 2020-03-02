/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_CARTPOLEDYNAMICS_H
#define MPCSIM_CARTPOLEDYNAMICS_H


#include "../../src/Dynamics.h"

/**
 * Implements the equations of motion of a Cart-Pole system.
 *
 * A classic system for studying controls is the cart-pole system, consisting of a cart on a linear track with pole
 * attached to the center with a frictionless joint. The pole has a large mass at the other end, and the control of the
 * system is in the horizontal force applied to the cart.
 *
 *                                                         O
 *                                                        /
 *                                                       /
 *                                              ________/____
 *                                             |       /     |
 *                                             |      o      |---> u
 *                                             |_____________|
 *                                    ____________O_______O____________
 *
 * The system is under-actuated, which makes the control challenging, especially in the presence of disturbances. A
 * problem of interest is "righting" the pole, which means bringing it from a "hanging down" configuration to a
 * "standing up" configuration, which is achieved by moving the cart back and forth quickly.
 */
class CartPoleDynamics : public Dynamics {

private:
    // Mass of the cart
    double m_mass_cart;

    // Mass of the pole
    double m_mass_pole;

    // Total mass of the cart and pole
    double m_mass_total;

    // Length of the pole;
    double m_length_pole;

public:
    /**
     * Constructor of the CartPoleDynamics class, that initialized the system parameters (masses, lengths, etc.)
     *
     * @param mass_cart
     * @param mass_pole
     * @param length_pole
     */
    CartPoleDynamics(double mass_cart, double mass_pole, double length_pole)
        : m_mass_cart(mass_cart),
          m_mass_pole(mass_pole),
          m_length_pole(length_pole)
    {
        // Compute the total mass of the cart and the pole
        this->m_mass_total = mass_cart + mass_pole;
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
     * Returns the number of elements (dimensionality) of the state vector, which is 4 for this system:
     *
     *  x = [z, z_dot, theta, theta_dot]^T
     *
     * @return 4
     */
    int getStateDimension() {return 4;}

    /**
     * Returns the number of elements (dimensionality) of the control vector, which is 1 for this system:
     *
     *  u = [u]
     *
     * @return 1
     */
    int getControlDimension() {return 1;}
};


#endif //MPCSIM_CARTPOLEDYNAMICS_H
