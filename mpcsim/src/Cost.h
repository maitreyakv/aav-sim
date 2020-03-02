/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_COST_H
#define MPCSIM_COST_H


#include <armadillo>

/**
 * Abstract class that implements the cost functions used in the computing of the optimal control
 *
 * This abstract class defines the terminal cost, the transition cost, and the gradients and hessians of these cost
 * functions. These are used by trajectory optimization algorithms in the Controller, like DDP, to determine the
 * optimal trajectory. "Optimal" in this sense, is entirely dependent on how these cost functions are defined, which
 * is why this is an abstract class. A user would need to define the cost functions by creating a class that inherits
 * this one and implements its methods.
 *
 * Note: These functions are entirely required by DDP, and may not be used by other trajectory optimizers. In the
 * future, if other optimizers are available, this class may be refactored to allow for integration with other methods
 * like PMP.
 */
class Cost {

public:
    /**
     * Terminal cost function, penalizes separation from the final state in a trajectory from the target state
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The terminal cost evaluated for the trajectory's final state
     */
    virtual double phi(arma::vec x_f, arma::vec x_star) = 0;

    /**
     * Gradient of the terminal cost function with respect to the state, evaluated at the terminal state
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The gradient of the terminal cost evaluated for the trajectory's final state
     */
    virtual arma::vec phi_x(arma::vec x_f, arma::vec x_star) = 0;

    /**
     * Hessian of the terminal cost function with respect to the state, evaluated at the terminal state
     *
     * @param x_f Vector which is the final state in a trajectory
     * @param x_star Vector which is the target state of the System
     * @return The hessian of the terminal cost evaluated for the trajectory's final state
     */
    virtual arma::mat phi_xx(arma::vec x_f, arma::vec x_star) = 0;

    /**
     * Transition cost of a state-action pair, multiplied by the algorithm time-step. Note, this time-step is time-step
     * used in discretizing time in the trajectory optimizer, such as DDP.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Transition cost evaulated for the state-action pair, multiplied by the time-step
     */
    virtual double L(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the state.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Gradient of the transition cost multiplied by the time-step, with respect to the state
     */
    virtual arma::vec L_x(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the control.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Gradient of the transition cost multiplied by the time-step, with respect to the control
     */
    virtual arma::vec L_u(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Hessian of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the state.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Hessian of the transition cost multiplied by the time-step, with respect to the state
     */
    virtual arma::mat L_xx(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Hessian of the transition cost of a state-action pair, multiplied by the algorithm time-step, with respect to
     * the control.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Hessian of the transition cost multiplied by the time-step, with respect to the control
     */
    virtual arma::mat L_uu(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Double gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with
     * respect to the state first and then the control.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Double gradient of the transition cost multiplied by the time-step, with respect to state, then control
     */
    virtual arma::mat L_xu(arma::vec x, arma::vec u, double dt) = 0;

    /**
     * Double gradient of the transition cost of a state-action pair, multiplied by the algorithm time-step, with
     * respect to the control first and then the state.
     *
     * @param x Vector containing the state to evaluate at.
     * @param u Vector containing the control to evaluate at.
     * @return Double gradient of the transition cost multiplied by the time-step, with respect to control, then state
     */
    virtual arma::mat L_ux(arma::vec x, arma::vec u, double dt) = 0;
};


#endif //MPCSIM_COST_H
