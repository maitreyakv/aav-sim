/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_CONTROLLER_H
#define MPCSIM_CONTROLLER_H


#include "Dynamics.h"
#include "Cost.h"

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

public:
    /**
     * Constructor for the Dynamics class that uses member-list initialization
     *
     * @param dynamics_ptr Pointer to the instance of a subclass of Dynamics that represents the dynamics of the system
     * @param cost_ptr Pointer to the instance of a subclass of Cost that represents the cost functions
     */
    Controller(Dynamics* dynamics_ptr, Cost* cost_ptr)
        : m_dynamics_ptr(dynamics_ptr),
          m_cost_ptr(cost_ptr)
    {}

    /**
     * Abstract method to compute the optimal control input of the System at the current state of the System. Whatever
     * method is used to determine this, which will implement a class that inherits from this one, will implement this
     * method using its own algorithm, i.e. DDP or PMP.
     * @return Vector with the optimal (or locally optimal) control input for the system at the current state
     */
    virtual arma::vec computeOptimalControl();
};


#endif //MPCSIM_CONTROLLER_H
