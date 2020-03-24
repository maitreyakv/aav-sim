/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "SamplerController.h"

arma::vec SamplerController::computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f) {
    // Compute the time-step of the time discretization
    double dt = (t_f - t_0) / (this->m_num_discretization - 1);

    // Allocate and compute the time-stamps of the discretized time horizon
    std::vector<double> t;
    for (int k = 0; k < this->m_num_discretization; k++) {
        t.push_back(t_0 + k * dt);
    }

    // Obtain state and control vector sizes
    int dim_x = this->m_dynamics_ptr->getStateDimension();
    int dim_u = this->m_dynamics_ptr->getControlDimension();

    // Allocate sequences for the trajectory
    std::vector<arma::vec> x(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));

    // Initialize the sequences with initial state
    x[0] = x_0;

    // Allocate sequences for the control sequence and the sampled change to the sequence
    std::vector<arma::vec> u(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));
    std::vector<arma::vec> u_best(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));
    std::vector<arma::vec> du(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));

    // Declare the best control sequences total cost and initialize it with the possible value
    double J_best = std::numeric_limits<double>::max();

    // Perform sampling of control sequences and simulate their trajectories
    for (int n = 0; n < this->m_num_samples; n++) {
        // Copy the best control sequence to the new sequence to be generated
        u = u_best;

        // Apply change to the best control sequence at every point in the sequence
        for (int k = 0; k < this->m_num_discretization; k++) {
            // Create random control change vector with elements in the range [0,1]
            for (int m = 0; m < dim_u; m++) {
                du[k](m) = rand() / (double) RAND_MAX;
            }

            // Scale the random vector for the change to 10% the range of the maximum control
            du[k] = 0.1 * (2.0 * this->m_u_max % du[k] - this->m_u_max);

            // Assign the control in the sequence
            u[k] = u[k] + du[k];
        }

        // Simulate trajectory using sampled control sequence
        for (int k = 0; k < this->m_num_discretization - 1; k++) {
            // Integrate state using Euler integration scheme
            x[k+1] = x[k] + this->m_dynamics_ptr->F(x[k], u[k], t[k]) * dt;
        }

        // Compute total cost of trajectory
        double J = this->m_cost_ptr->phi(x.back(), x_star);
        for (int k = 0; k < this->m_num_discretization; k++) {
            J = J + this->m_cost_ptr->L(x[k], u[k], dt);
        }

        // Replace the best control sequence with the current sequence if the cost is lower than the best cost
        if (J < J_best && !std::isnan(J)) {
            // Update the best control sequence
            u_best = u;

            // Update the best cost
            J_best = J;
        }
    }

    // Return the first control input in the best control sequence
    return u_best[0];
}