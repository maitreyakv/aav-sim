/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_SAMPLERCONTROLLER_H
#define MPCSIM_SAMPLERCONTROLLER_H


#include <armadillo>

#include "Controller.h"


class SamplerController : public Controller {

private:
    // Number of time discretizations used in the sampled trajectory
    int m_num_discretization;

    // Number of samples to generate
    int m_num_samples;

public:
    /**
     * Constructor for the SamplerController class that uses member-list initialization
     *
     * @param dynamics_ptr Pointer to the Dynamics object used by the System using this Controller
     * @param cost_ptr Pointer to the Cost object used by this Controller
     * @param u_max Vector with the maximum control magnitudes
     * @param num_discretization Number of discretizations of the time horizon in sampled trajectories
     * @param num_samples Number of samples to generate
     */
    SamplerController(Dynamics* dynamics_ptr, Cost* cost_ptr, arma::vec u_max, int num_discretization, int num_samples)
        : Controller(dynamics_ptr, cost_ptr, u_max),
          m_num_discretization(num_discretization),
          m_num_samples(num_samples)
    {}


    arma::vec computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f);
};


#endif //MPCSIM_SAMPLERCONTROLLER_H
