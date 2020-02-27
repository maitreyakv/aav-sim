/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "DDPController.h"

arma::vec DDPController::computeOptimalControl(arma::vec x_0, double t_0, double t_f) {
    // Compute the time-step of the time discretization used in the algorithm
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
    std::vector<arma::vec> x_new(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));

    // Initialize the sequences with
    x[0] = x_0;
    x_new[0] = x_0;

    // Allocate sequence for the control
    std::vector<arma::vec> u;

    // Initialize control with zeros
    for (int k = 0; k < this->m_num_discretization; k++) {
        u.push_back(arma::zeros<arma::vec>(dim_u));
    }

    // Allocate sequences for value function and its gradient and Hessian
    std::vector<double> V(this->m_num_discretization, 0.0);
    std::vector<arma::vec> V_x(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::mat> V_xx(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Allocate sequences for the state-action value function's various first and second derivatives
    std::vector<arma::vec> Q_x(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::vec> Q_u(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::mat> Q_xx(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_uu(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_xu(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_ux(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Generate initial trajectory using initial random control sequence
    for (int k = 0; k < this->m_num_discretization - 1; k++) {
        // Integrate state using Euler integration scheme
        x_new[k+1] = x_new[k] + this->m_dynamics_ptr->F(x_new[k], u[k], t[k]) * dt;
    }

    // TODO: Perform DDP to Optimize u

    // Return the first control input in the optimal control sequence
    return u[0];
}