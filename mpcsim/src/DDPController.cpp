/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "DDPController.h"

arma::vec DDPController::computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f) {
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

    // Initialize control with random control
    for (int k = 0; k < this->m_num_discretization; k++) {
        u.push_back({rand() / (double) RAND_MAX,
                     rand() / (double) RAND_MAX,
                     rand() / (double) RAND_MAX,
                     rand() / (double) RAND_MAX});
        u[k] = 2.0 * this->m_u_max % u[k] - this->m_u_max;
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

    // Allocate vectors for the feed-back and feed-forward terms in the control update
    arma::vec du_ff = arma::zeros<arma::vec>(dim_u);
    arma::vec du_fb = arma::zeros<arma::vec>(dim_u);

    // Perform DDP algorithm by performing forward and backwards passes
    for (int i = 0; i < this->m_num_iteration; i++) {
        // Update the control if not the first iteration
        if (i > 0) {
            // Compute feed-forward and feed-backward terms and update control (Forward Pass)
            for (int k = 0; k < this->m_num_discretization - 1; k++) {
                // Compute feed-forward and feed-backward terms
                du_ff = - arma::solve(Q_uu[k], Q_u[k]);
                du_fb = - arma::solve(Q_uu[k], Q_ux[k]) * (x_new[k] - x[k]);

                // TODO: Implement control clamping

                // Update control using correction terms scaled by the learning rate
                u[k] = u[k] + this->m_learning_rate * (du_ff + du_fb);

                // Propagate the trajectory using the updated control
                x_new[k+1] = x_new[k] + this->m_dynamics_ptr->F(x_new[k], u[k], t[k]) * dt;
            }
        }

        // Copy the new trajectory into the old trajectory for the purpose of comparing trajectories between iterations
        x = x_new;

        // Compute total cost of trajectory
        double J = this->m_cost_ptr->phi(x.back(), x_star);
        for (int k = 0; k < this->m_num_discretization; k++) {
            J = J + this->m_cost_ptr->L(x[k], u[k], dt);
        }

        // Compute the value function, its gradient and Hessian of the last state using the terminal cost
        V.back() = this->m_cost_ptr->phi(x.back(), x_star);
        V_x.back() = this->m_cost_ptr->phi_x(x.back(), x_star);
        V_xx.back() = this->m_cost_ptr->phi_xx(x.back(), x_star);

        // Propagate the state-action value function and its derivatives backwards (Backwards Pass)
        for (int k = this->m_num_discretization - 2; k >=0 ; k--) {
            // Compute state-action value function and its derivatives
            Q_x[k] = this->m_cost_ptr->L_x(x[k], u[k], dt) + this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_x[k+1];
            Q_u[k] = this->m_cost_ptr->L_u(x[k], u[k], dt) + this->m_dynamics_ptr->Beta(x[k], u[k], dt).t() * V_x[k+1];
            Q_xx[k] = this->m_cost_ptr->L_xx(x[k], u[k], dt)
                    + this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_xx[k+1]
                    * this->m_dynamics_ptr->Phi(x[k], u[k], dt);
            Q_uu[k] = this->m_cost_ptr->L_uu(x[k], u[k], dt)
                    + this->m_dynamics_ptr->Beta(x[k], u[k], dt).t() * V_xx[k+1]
                    * this->m_dynamics_ptr->Beta(x[k], u[k], dt);
            Q_xu[k] = this->m_cost_ptr->L_xu(x[k], u[k], dt)
                    + this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_xx[k+1]
                    * this->m_dynamics_ptr->Beta(x[k], u[k], dt);
            Q_ux[k] = this->m_cost_ptr->L_ux(x[k], u[k], dt)
                    + this->m_dynamics_ptr->Beta(x[k], u[k], dt).t()
                    * V_xx[k+1] * this->m_dynamics_ptr->Phi(x[k], u[k], dt);

            // Compute the value function gradient and Hessian during the backwards pass
            V_x[k] = Q_x[k] - Q_xu[k] * arma::solve(Q_uu[k], Q_u[k]);
            V_xx[k] = Q_xx[k] - Q_xu[k] * arma::solve(Q_uu[k], Q_ux[k]);
        }
    }

    // Return the first control input in the optimal control sequence
    return u[0];
}