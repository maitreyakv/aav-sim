/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_SIMULATION_H
#define MPCSIM_SIMULATION_H


#include "System.h"

/**
 * Implementation of a single simulation in the program
 *
 * This class represents a single simulation, which is used to organize the different elements of the simulation, as
 * well as provide a single object in which to initialize, execute, and export results from a simulation.
 */
class Simulation {

private:
    // System to be simulated
    System m_system;

public:
    /**
     * Constructor for a Simulation.
     */
    Simulation();

    /**
     * Executes the simulation process by integrating the dynamics of the system while computing and using the optimal
     * control input
     *
     * @param time Total time to simulate the system for from the initial state
     */
    void simulate(double time);
};


#endif //MPCSIM_SIMULATION_H
