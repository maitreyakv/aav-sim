/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAVSIM_QUADRATICOBSTACLECOST_H
#define AAVSIM_QUADRATICOBSTACLECOST_H


#include <armadillo>

#include "Cost.h"

/**
 * Implementation of quadratic cost function with obstacle avoidance.
 *
 * The quadratic cost function for the control and terminal cost are the same as in the QuadraticCost class. This
 * class also implemented Gaussian-form cost functions to penalize the proximity of the state to a set of obstacle
 * states. This in effect, drives the optimized trajectories away from these states, which in the case of a physical
 * obstacle, means distancing a vehicle from a set of obstacles.
 *
 * The transition cost is multiplied by the algorithm time-step, since the cost is being applied during the entire
 * transition of the state given the state-action pair (x,u).
 */
class QuadraticObstacleCost : public Cost {

private:
    // Matrix of weights used in quadratic terminal cost
    arma::mat m_Q_f;

    // Matrix of weights used in the quadratic transition cost, with only control cost
    arma::mat m_R;

    // List of obstacles
    std::vector<arma::vec> m_obstacles;

    // Weight matrix penalizing distance to obstacles
    arma::mat m_sigma;

public:
    /**
     * Constructor using member-list initialization
     *
     * @param Q_f Matrix of weights used in quadratic terminal cost
     * @param R Matrix of weights used in the quadratic transition cost, with only control cos
     * @param obstacles List of obstacles state vectors
     * @param sigma Penalizing matrix for the obstacles
     */
    QuadraticObstacleCost(arma::mat Q_f, arma::mat R, std::vector<arma::vec> obstacles, arma::mat sigma)
        : m_Q_f(Q_f),
          m_R(R),
          m_obstacles(obstacles),
          m_sigma(sigma)
    {}

    /**
     * Terminal cost function of the form
     *
     *      phi = 0.5 * (x_f - x_star)^T * Q_f * (x_f - x_star)
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The terminal cost evaluated for the trajectory's final state
     */
    double phi(arma::vec x_f, arma::vec x_star);

    /**
     * Gradient of the terminal cost function with respect to the state, of the form
     *
     *      phi_x = Q_f * (x_f - x_star)
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The gradient of the terminal cost evaluated for the trajectory's final state
     */
    arma::vec phi_x(arma::vec x_f, arma::vec x_star);

    /**
     * Hessian of the terminal cost function with respect to the state, of the form
     *
     *      phi_xx = Q_f
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The hessian of the terminal cost evaluated for the trajectory's final state
     */
    arma::mat phi_xx(arma::vec x_f, arma::vec x_star);

    /**
     * Transition cost of a state-action pair, multiplied by the algorithm time-step, of the form
     *
     *      L = 0.5 * u^T * R * u * dt
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Transition cost evaulated for the state-action pair, multiplied by the time-step
     */
    double L(arma::vec x, arma::vec u, double dt);

    /**
     * Gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the state, of the form
     *
     *      L_x = [0]_Nx1
     *
     * where N is the size of the state vector
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Gradient of the transition cost multiplied by the time-step, with respect to the state
     */
    arma::vec L_x(arma::vec x, arma::vec u, double dt);

    /**
     * Gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the control, of the form
     *
     *      L_u = R * u * dt
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Gradient of the transition cost multiplied by the time-step, with respect to the control
     */
    arma::vec L_u(arma::vec x, arma::vec u, double dt);

    /**
     * Hessian of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the state, of the form
     *
     *      L_xx = [0]_NxN
     *
     * where N is the size of the state vector
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Hessian of the transition cost multiplied by the time-step, with respect to the state
     */
    arma::mat L_xx(arma::vec x, arma::vec u, double dt);

    /**
     * Hessian of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the control, of the form
     *
     *      L_uu = R * dt
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Hessian of the transition cost multiplied by the time-step, with respect to the control
     */
    arma::mat L_uu(arma::vec x, arma::vec u, double dt);

    /**
     * Double gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with
     * respect to the state first and then the control, of the form
     *
     *      L_xu = [0]_NxM
     *
     * where N is the size of the state vector, and M is the size of the control vector
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Double gradient of the transition cost multiplied by the time-step, with respect to state, then control
     */
    arma::mat L_xu(arma::vec x, arma::vec u, double dt);

    /**
     * Double gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with
     * respect to the control first and then the state, of the form
     *
     *      L_ux = [0]_MxN
     *
     * where N is the size of the state vector, and M is the size of the control vector
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Double gradient of the transition cost multiplied by the time-step, with respect to control, then state
     */
    arma::mat L_ux(arma::vec x, arma::vec u, double dt);
};


#endif //AAVSIM_QUADRATICOBSTACLECOST_H
