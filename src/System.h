/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_SYSTEM_H
#define MPCSIM_SYSTEM_H


#include <armadillo>

#include "Controller.h"
#include "Dynamics.h"

/**
 * Implementation of a System in the Simulation
 *
 * The words System in this context refers to the collection of:
 *      1. The Dynamics of the System (equations of motion, governing equations, etc.)
 *      2. The Controller that computed the optimal control (DDP, PMP, etc.) with some frequency
 * The System therefore collects the physical process being controlled (Dynamics), as well as the algorithm that
 * computes the optimal control (Controller) in one object. The object is also used by "odeint" to produce the state
 * derivatives as a function of some state.
 */
class System {

private:
    // Dynamics of the System
    Dynamics* m_Dynamics;

    // Controller of the System
    Controller* m_Controller;

    // Current control input being used
    arma::vec u;

public:
    /**
     * Overload the () operator using the equations of motion for the purpose of using "odeint" to perform the
     * integration of the equations of motion. The system of the simulation is of the form
     *
     *      dx/dt = F(x,u)
     *
     * which means that the () operator is equivalent to F(x,u), where "u" is the current control input, which is a
     * member of this class, and is updated with some frequency by the Simulation. In other words, it is the
     * function that evaluates the state derivative using a state of the system.
     *
     * @param x Address of a vector containing a state where the state derivative is to be evaluated
     * @param dxdt Address of a vector where the calculated state derivative will be placed in memory
     * @param t The time at which the state derivative will be evaluated, if the dynamics are explicit in time
     */
    void operator()  (const arma::vec& x , arma::vec& dxdt , const double t);

    /**
     * Uses the Controller to recompute the best control input given the current state of the system. This process
     * is repeated with some frequency during the simulation, which what the MPC method does, i.e. continually recompute
     * the optimal trajectory with a moving time horizon as the system evolves, in order to account for disturbances.
     */
    void updateControl();
};


#endif //MPCSIM_SYSTEM_H
