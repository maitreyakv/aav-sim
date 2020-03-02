/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_DDPCONTROLLER_H
#define MPCSIM_DDPCONTROLLER_H


#include <armadillo>

#include "Controller.h"

/**
 * Implements the DDP Algorithm as a Controller for MPC
 *
 * The Differential Dynamic Programming algorithm is a powerful method in optimal control that computes a locally
 * optimal control sequence for a finite time horizon. Because Dynamic Programming is infeasible for continuous
 * systems, due to the "curse of dimensionality", Differential Dynamic Programming performs Dynamic Programming
 * locally around some nominal trajectory and corresponding control sequence, in order to adjust the control sequence
 * to a locally optimal one. The algorithm performs a series of "forward" and "backward" passes on a nominal trajectory,
 * adjusting the control at each of the discretized states in the time horizon to reduce some user-defined cost
 * function.
 *
 * Note: The time-step used in the DDP algorithm is NOT the same as the time-step used in the Simulation. The step used
 * in DDP is much larger, since the accuracy of the computed locally optimal trajectory is not important, due to the
 * fact that in MPC, only the first control input in the computed sequence is implemented. The step in the actual
 * simulation is much larger, since we are interested in determining in an accurate manner how the system would evolve,
 * so we use a small step to reduce the introduction of truncation errors into the system.
 */
class DDPController : public Controller {

private:
    // Number of time discretizations used in DDP
    int m_num_discretization;

    // Number of iterations in forward-backward passes in DDP
    int m_num_iteration;

    // Learning rate in control update in DDP
    double m_learning_rate;

public:
    /**
     * Constructor for the DDPController class that uses member-list initialization
     *
     * @param dynamics_ptr Pointer to the Dynamics object used by the System using this Controller
     * @param cost_ptr Pointer to the Cost object used by this Controller
     * @param u_max Vector with the maximum control magnitudes
     * @param num_discretization Number of discretizations of the time horizon in the DDP algorithm
     * @param num_iteration Number of iterations (forward-backward-passes) used by the algorithm
     * @param learning_rate Learning rate to use when updating the control each iteration
     */
    DDPController(Dynamics* dynamics_ptr, Cost* cost_ptr, arma::vec u_max, int num_discretization, int num_iteration,
                  double learning_rate)
        : Controller(dynamics_ptr, cost_ptr, u_max),
          m_num_discretization(num_discretization),
          m_num_iteration(num_iteration),
          m_learning_rate(learning_rate)
    {}

    /**
     * Uses the DDP algorithm to compute the optimal control sequence for the finite time horizon, and the then returns
     * the first control element as the selected optimal control for the current MPC step. The DDP method involves:
     *
     *      1. Initializing control sequence randomly, or with zero control
     *      2. Propagating trajectory using control sequence (forward pass)
     *      3. Evaluating value function, its gradient and Hessian, at the final state of the trajectory
     *      4. Propagating these three values backwards to the initial state (backward pass)
     *      5. Update the control with a computed correction
     *      6. Repeat from step 2 until convergence on locally optimal trajectory
     *
     * @param x_0 Initial state of the time horizon, which is the current state of the System
     * @param t_0 Initial time of the time horizon, which is the current time of the System
     * @param x_star Target state to attempt to reach in the horizon
     * @param t_f Final time of the time horizon, which is t_0 + length of the time horizon
     * @return Vector containing the first element in the (locally) optimal control sequence
     */
    arma::vec computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f);
};


#endif //MPCSIM_DDPCONTROLLER_H
