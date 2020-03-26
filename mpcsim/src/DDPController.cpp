/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "DDPController.h"

arma::vec DDPController::computeOptimalControl(arma::vec x_0, double t_0, arma::vec x_star, double t_f) {
    // Set up matrix system solver options for armadillo
    arma::solve_opts::opts arma_solver_opts = arma::solve_opts::likely_sympd
                                            + arma::solve_opts::refine
                                            + arma::solve_opts::allow_ugly;

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

    // Initialize the sequences with initial state
    x[0] = x_0;
    x_new[0] = x_0;

    // Allocate sequence for the nominal control sequence and the updated control sequence
    std::vector<arma::vec> u(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));
    std::vector<arma::vec> u_new(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));

    // Initialize control with random control no more than 10% of maximum allowable control
    for (int k = 0; k < this->m_num_discretization; k++) {
        // Create random control vector with elements in the range [0,1]
        arma::vec u_rand = arma::zeros<arma::vec>(dim_u);
        for (int m = 0; m < dim_u; m++) {
            u_rand(m) = rand() / (double) RAND_MAX;
        }

        // Scale the random vector to 10% of maximum control
        u[k] = 0.1 * (2.0 * this->m_u_max % u_rand - this->m_u_max);
    }

    // Allocate sequences for value function and its gradient and Hessian
    //std::vector<double> V(this->m_num_discretization, 0.0);
    std::vector<arma::vec> V_x(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::mat> V_xx(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Allocate sequence for the expected changes in the value function
    std::vector<double> delta_V(this->m_num_discretization, 0.0);

    // Allocate sequences for the state-action value function's various first and second derivatives
    std::vector<arma::vec> Q_x(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::vec> Q_u(this->m_num_discretization, arma::zeros<arma::vec>(dim_x));
    std::vector<arma::mat> Q_xx(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_uu(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_xu(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_ux(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Allocate sequences for the feed-forward and feed-back gains
    std::vector<arma::vec> gain_ff(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));
    std::vector<arma::mat> gain_fb(this->m_num_discretization, arma::zeros<arma::mat>(dim_u, dim_x));

    // Allocate the regularization parameter
    double mu;

    // Generate initial trajectory using initial random control sequence
    for (int k = 0; k < this->m_num_discretization - 1; k++) {
        // Integrate state using Euler integration scheme
        x[k+1] = x[k] + this->m_dynamics_ptr->F(x[k], u[k], t[k]) * dt;

#ifdef DEBUG
        for (int i = 0; i < dim_x; i++) {
            if (std::isnan(x[k+1](i))) {
                std::cout << "error: DDP initialized trajectory contains NaN\n"; exit(-1);
            }
        }
#endif
    }

    // Allocate vectors for the feed-back and feed-forward terms in the control update
    arma::vec du_ff = arma::zeros<arma::vec>(dim_u);
    arma::vec du_fb = arma::zeros<arma::vec>(dim_u);

    // Perform DDP algorithm by performing forward and backwards passes
    for (int i = 0; i < this->m_num_iteration; i++) {
        // Compute total cost of current trajectory
        double J = this->m_cost_ptr->phi(x.back(), x_star);
        for (int k = 0; k < this->m_num_discretization; k++) {
            J = J + this->m_cost_ptr->L(x[k], u[k], dt);
        }

        // Compute the value function, its gradient and Hessian of the last state using the terminal cost
        //V.back() = this->m_cost_ptr->phi(x.back(), x_star);
        V_x.back() = this->m_cost_ptr->phi_x(x.back(), x_star);
        V_xx.back() = this->m_cost_ptr->phi_xx(x.back(), x_star);

        // Perform the backward pass
        for (int k = this->m_num_discretization - 2; k >= 0; k--) {
            // Compute state-action value function and its derivatives
            Q_x[k] = this->m_cost_ptr->L_x(x[k], u[k], dt) +
                     this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_x[k + 1];
            Q_u[k] = this->m_cost_ptr->L_u(x[k], u[k], dt) +
                     this->m_dynamics_ptr->Beta(x[k], u[k], dt).t() * V_x[k + 1];
            Q_xx[k] = this->m_cost_ptr->L_xx(x[k], u[k], dt)
                      + this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_xx[k + 1]
                        * this->m_dynamics_ptr->Phi(x[k], u[k], dt);
            Q_uu[k] = this->m_cost_ptr->L_uu(x[k], u[k], dt)
                      + this->m_dynamics_ptr->Beta(x[k], u[k], dt).t() * V_xx[k + 1]
                        * this->m_dynamics_ptr->Beta(x[k], u[k], dt);
            Q_xu[k] = this->m_cost_ptr->L_xu(x[k], u[k], dt)
                      + this->m_dynamics_ptr->Phi(x[k], u[k], dt).t() * V_xx[k + 1]
                        * this->m_dynamics_ptr->Beta(x[k], u[k], dt);
            Q_ux[k] = this->m_cost_ptr->L_ux(x[k], u[k], dt)
                      + this->m_dynamics_ptr->Beta(x[k], u[k], dt).t()
                        * V_xx[k + 1] * this->m_dynamics_ptr->Phi(x[k], u[k], dt);

            // Ensure that Q_uu is numerically exactly symmetric
            Q_uu[k] = 0.5 * ( arma::symmatu(Q_uu[k]) + arma::symmatl(Q_uu[k]) );

            // Obtain eigenvalues of unregularized Q_uu matrix
            arma::vec Q_uu_eigval = arma::eig_sym(Q_uu[k]);

            // Compute regularization parameter using smallest eigenvalue if Q_uu is not positive definite
            mu = (Q_uu_eigval[0] < 0.0) ? -Q_uu_eigval[0] : 0.0;

            // Regularize state-action value functions
            Q_uu[k] = Q_uu[k] + mu * arma::eye(dim_u, dim_u);

            // Compute the feed-forward and feed-backward gains
            gain_ff[k] = -arma::solve(Q_uu[k], Q_u[k], arma_solver_opts);
            gain_fb[k] = -arma::solve(Q_uu[k], Q_ux[k], arma_solver_opts);

            // Compute the value function gradient and Hessian during the backwards pass
            V_x[k] = Q_x[k] - Q_xu[k] * arma::solve(Q_uu[k], Q_u[k], arma_solver_opts);
            V_xx[k] = Q_xx[k] - Q_xu[k] * arma::solve(Q_uu[k], Q_ux[k], arma_solver_opts);

            // Compute the expected change in the value function
            delta_V[k] = arma::as_scalar( -0.5 * Q_u[k].t() * arma::solve(Q_uu[k], Q_u[k], arma_solver_opts) );
        }

        // Compute feed-forward and feed-backward terms and update control (Forward Pass)
        for (int k = 0; k < this->m_num_discretization - 1; k++) {
            // Compute feed-forward and feed-backward terms
            du_ff = gain_ff[k];
            du_fb = gain_fb[k] * (x_new[k] - x[k]);

            // Limit feed-forward component of control update using simple clamping
            for (int m = 0; m < dim_u; m++) {
                du_ff(m) = fmin(this->m_u_max(m), fmax(-this->m_u_max(m), du_ff(m) + u[k](m))) - u[k](m);
            }

            // Update control using correction terms scaled by the learning rate
            u_new[k] = u[k] + this->m_learning_rate * du_ff + du_fb;

            // Propagate the trajectory using the updated control
            x_new[k + 1] = x_new[k] + this->m_dynamics_ptr->F(x_new[k], u_new[k], t[k]) * dt;
        }

        // Compute total cost of new trajectory
        double J_new = this->m_cost_ptr->phi(x_new.back(), x_star);
        for (int k = 0; k < this->m_num_discretization; k++) {
            J_new = J_new + this->m_cost_ptr->L(x_new[k], u_new[k], dt);
        }

        // Copy the updated control into the current nominal control
        u = u_new;

        // Copy the new trajectory into the old trajectory for the purpose of comparing trajectories between iterations
        x = x_new;

        // TEMP
        std::cout << J << "," << J_new << std::endl;
    }

    // Return the first control input in the optimal control sequence
    return u[0];
}