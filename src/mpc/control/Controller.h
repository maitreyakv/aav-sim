/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAVSIM_CONTROLLER_H
#define AAVSIM_CONTROLLER_H


#include "../../simulation/Dynamics.h"
#include "../cost/Cost.h"

/**
 * Implements the controlling algorithm for the System.
 *
 * This class implements a control algorithm, like DDP, PMP, or some other method, that is called by the MPC algorithm
 * at each of its steps in order to "continually" update the control to account for disturbances. The Controller takes
 * the current state of the system, and finds the optimal, or locally optimal, control sequence to generate a
 * trajectory to take the state as close as possible within a finite time horizon. This horizon keeps moving each time
 * this process repeats (receding horizon control).
 */
class Controller {

protected:
    // Pointer to Dynamics of the System being controlled
    Dynamics* m_dynamics_ptr;

    // Pointer to Cost functions that were defined to tune the Controller
    Cost* m_cost_ptr;

    // Maximum control magnitude
    arma::vec m_u_max;

public:
    /**
     * Constructor for the Controller class that uses member-list initialization
     *
     * @param dynamics_ptr Pointer to the instance of a subclass of Dynamics that represents the dynamics of the system
     * @param cost_ptr Pointer to the instance of a subclass of Cost that represents the cost functions
     * @param u_max Vector with the maximum magnitude of the control inputs
     */
    Controller(Dynamics* dynamics_ptr, Cost* cost_ptr, arma::vec u_max)
        : m_dynamics_ptr(dynamics_ptr),
          m_cost_ptr(cost_ptr),
          m_u_max(u_max)
    {}

    /**
     * Abstract method to compute the optimal control input of the System at the current state of the System. Whatever
     * method is used to determine this, which will implement a class that inherits from this one, will implement this
     * method using its own algorithm, i.e. DDP or PMP.
     *
     * @param x_0 Initial state of the time horizon, which is the current state of the System
     * @param t_0 Initial time of the time horizon, which is the current time of the System
     * @param x_star Target state to attempt to reach in the horizon
     * @param t_f Final time of the time horizon, which is t_0 + length of the time horizon
     * @return Vector containing the optimal control to implement at this step in MPC
     */
    virtual arma::vec computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f) = 0;
};


#endif //AAVSIM_CONTROLLER_H
