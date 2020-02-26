/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_DYNAMICS_H
#define MPCSIM_DYNAMICS_H


#include <armadillo>

/**
 * Abstract implementation of the Dynamics of a System.
 *
 * This abstract class is used to represent the equations of motion of a system, in state-space form. For a user to be
 * able to use this program to simulate any system of their choice, they would need to create a class that inherits
 * from this one, and that implements the equations of motion in state-space form
 *
 *      dx/dt = F(x,u)
 *
 * for their particular system. In addition, in order to use DDP as the control-determining algorithm in each iteration
 * of the MPC algorithm, the user will also need to implement the two gradients of the dynamics as well, so that the
 * DDP algorithm has access to the linearized dynamics.
 *
 * Tip: If the closed form equations of motion for "F(x,u)" are known, then these gradients with respect to "x" and
 * "u" can be easily computed using a linear algebra software like "SciPy" in Python
 */
class Dynamics {

public:
    /**
     * Abstract method to compute the state derivative dxdt = F(x,u) of the Dynamics
     *
     * @param x Vector with the state "x" to evaluate F(x,u) at
     * @param u Vector with the control input "u" to evaluate F(x,u) at
     * @param t The time to evaluate F(x,u) at, if it is explicit in time, i.e. if dx/dt = F(x,u,t)
     * @return Vector containing the state derivative dx/dt = F(x,u)
     */
    virtual arma::vec F(arma::vec x, arma::vec u, double t) = 0;

    /**
     * Abstract method to compute the gradient of the state derivative F_x of the Dynamics with respect to the state,
     * returned in the form
     *
     *      Phi = I + F_x(x,u) * dt
     *
     * This is specifically required by Differential Dynamic Programming, as it is used to linearize the dynamics of
     * the system in the algorithm.
     *
     * @param x Vector with the state "x" to evaluate F_x(x,u) at
     * @param u Vector with the control input "u" to evaluate F_x(x,u) at
     * @param dt The discretized time step in the DDP algorithm, t(k+1) - t(k)
     * @return Vector containing Phi
     */
    virtual arma::vec Phi(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Abstract method to compute the gradient of the state derivative F_x of the Dynamics with respect to the control,
     * returned in the form
     *
     *      Beta = F_u(x,u) * dt
     *
     * This is specifically required by Differential Dynamic Programming, as it is used to linearize the dynamics of
     * the system in the algorithm.
     *
     * @param x Vector with the state "x" to evaluate F_x(x,u) at
     * @param u Vector with the control input "u" to evaluate F_x(x,u) at
     * @param dt The discretized time step in the DDP algorithm, t(k+1) - t(k)
     * @return Vector containing Beta
     */
    virtual arma::vec Beta(arma::vec x, arma::vec u, double dt) = 0;
};


#endif //MPCSIM_DYNAMICS_H
