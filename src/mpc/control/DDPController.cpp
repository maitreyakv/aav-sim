/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "DDPController.h"

// TODO: add proper comments/documentation and clean code


arma::mat contract_vector_and_tensor(arma::vec vect, arma::cube tens) {
    arma::mat result = arma::zeros<arma::mat>(tens.n_cols, tens.n_slices);
    for (int i = 0; i < (int) tens.n_cols; i++) {
        for (int j = 0; j < (int) tens.n_slices; j++) {
            for (int k = 0; k < (int) vect.n_elem; k++) {
                result(i,j) = result(i,j) + vect(k) * tens(k,i,j);
            }
        }
    }
    return result;
}


// TEMP:
const double delta_0 = 2.0;
const double mu_min = pow(10, -6);
const int max_num_backward_pass = 10;
const int max_num_forward_pass = 10;

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

    if (this->m_use_initial_control) {
        u = this->m_u_initial;
    }
    else {
        // Initialize control with random control no more than 10% of maximum allowable control
        for (int k = 0; k < this->m_num_discretization; k++) {
            // Create random control vector with elements in the range [0,1]
            arma::vec u_rand = arma::zeros<arma::vec>(dim_u);
            for (int m = 0; m < dim_u; m++) {
                // TEMP
                u_rand(m) = 0.0 * rand() / (double) RAND_MAX;
            }

            // Scale the random vector to 10% of maximum control
            u[k] = 0.1 * (2.0 * this->m_u_max % u_rand - this->m_u_max);
        }

        this->m_use_initial_control = true;
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
    //std::vector<arma::mat> Q_xu(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_ux(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Allocate sequences for the regularized state-action value function derivatives
    std::vector<arma::mat> Q_uu_reg(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));
    std::vector<arma::mat> Q_ux_reg(this->m_num_discretization, arma::zeros<arma::mat>(dim_x, dim_x));

    // Allocate sequences for the feed-forward and feed-back gains
    std::vector<arma::vec> gain_ff(this->m_num_discretization, arma::zeros<arma::vec>(dim_u));
    std::vector<arma::mat> gain_fb(this->m_num_discretization, arma::zeros<arma::mat>(dim_u, dim_x));

    // Allocate the regularization parameter and modification factor
    double mu = 0.0;
    double delta = 1.0;

    // Generate initial trajectory using initial random control sequence
    for (int k = 0; k < this->m_num_discretization - 1; k++) {
        // Integrate state using Euler integration scheme
        x[k+1] = this->m_dynamics_ptr->f(x[k], u[k], t[k], dt);

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
        if (i > 0) {
            for (int k = 0; k < this->m_num_discretization; k++) {
                J = J + this->m_cost_ptr->L(x[k], u[k], dt);
            }
        }

        // Compute the value function, its gradient and Hessian of the last state using the terminal cost
        //V.back() = this->m_cost_ptr->phi(x.back(), x_star);
        V_x.back() = this->m_cost_ptr->phi_x(x.back(), x_star);
        V_xx.back() = this->m_cost_ptr->phi_xx(x.back(), x_star);

        // Initialize counter for number of backward passes
        int num_backward_pass = 0;

        // Initialize a boolean for if all the Q_uu matrices are positive definite
        bool is_Q_uu_reg_positive_definite;

        // Perform the backward pass repeatedly while using regularization until Q_uu is positive definite
        do {
            // Initialize the boolean for positive definite check to be true
            is_Q_uu_reg_positive_definite = true;

            // Perform backward pass
            for (int k = this->m_num_discretization - 2; k >= 0; k--) {
                // Compute state-action value function and its derivatives
                Q_x[k] = this->m_cost_ptr->L_x(x[k], u[k], dt) +
                         this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt).t() * V_x[k + 1];
                Q_u[k] = this->m_cost_ptr->L_u(x[k], u[k], dt) +
                         this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt).t() * V_x[k + 1];
                Q_xx[k] = this->m_cost_ptr->L_xx(x[k], u[k], dt)
                        + this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt).t() * V_xx[k + 1]
                        * this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt)
                        + contract_vector_and_tensor(V_x[k+1], this->m_dynamics_ptr->f_xx(x[k], u[k], t[k], dt));
                Q_uu[k] = this->m_cost_ptr->L_uu(x[k], u[k], dt)
                        + this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt).t() * V_xx[k + 1]
                        * this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt)
                        + contract_vector_and_tensor(V_x[k+1], this->m_dynamics_ptr->f_uu(x[k], u[k], t[k], dt));
                //Q_xu[k] = this->m_cost_ptr->L_xu(x[k], u[k], dt)
                //          + this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt).t() * V_xx[k + 1]
                //            * this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt);
                Q_ux[k] = this->m_cost_ptr->L_ux(x[k], u[k], dt)
                        + this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt).t()
                        * V_xx[k + 1] * this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt)
                        + contract_vector_and_tensor(V_x[k+1], this->m_dynamics_ptr->f_ux(x[k], u[k], t[k], dt));

                // Compute the regularized state-action value function derivatives
                Q_uu_reg[k] = this->m_cost_ptr->L_uu(x[k], u[k], dt)
                            + this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt).t()
                            * (V_xx[k + 1] + mu * arma::eye(dim_x, dim_x))
                            * this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt)
                            + contract_vector_and_tensor(V_x[k+1], this->m_dynamics_ptr->f_uu(x[k], u[k], t[k], dt));
                Q_ux_reg[k] = this->m_cost_ptr->L_ux(x[k], u[k], dt)
                            + this->m_dynamics_ptr->f_u(x[k], u[k], t[k], dt).t()
                            * (V_xx[k + 1] + mu * arma::eye(dim_x, dim_x))
                            * this->m_dynamics_ptr->f_x(x[k], u[k], t[k], dt)
                            + contract_vector_and_tensor(V_x[k+1], this->m_dynamics_ptr->f_ux(x[k], u[k], t[k], dt));

                // Ensure that the regularied Q_uu is numerically symmetric
                Q_uu_reg[k] = 0.5 * (arma::symmatl(Q_uu_reg[k]) + arma::symmatu(Q_uu_reg[k]));

                // Check if the regularized Q_uu is positive definite
                if (!Q_uu_reg[k].is_sympd()) {
                    // Set the flag to false to trigger another backward pass attempt
                    is_Q_uu_reg_positive_definite = false;
                }

                // Compute the feed-forward and feed-backward gains
                gain_ff[k] = -arma::solve(Q_uu_reg[k], Q_u[k], arma_solver_opts);
                gain_fb[k] = -arma::solve(Q_uu_reg[k], Q_ux_reg[k], arma_solver_opts);

                // Compute the value function gradient and Hessian during the backwards pass
                V_x[k] = Q_x[k] + gain_fb[k].t() * Q_uu[k] * gain_ff[k] + gain_fb[k].t() * Q_u[k]
                         + Q_ux[k].t() * gain_ff[k];
                V_xx[k] = Q_xx[k] + gain_fb[k].t() * Q_uu[k] * gain_fb[k] + +gain_fb[k].t() * Q_ux[k]
                          + Q_ux[k].t() * gain_fb[k];

                // Compute the expected change in the value function
                delta_V[k] = arma::as_scalar(0.5 * gain_ff[k].t() * Q_uu[k].t() * gain_ff[k] + gain_ff[k].t() * Q_u[k]);
            }

            // Decrease the regularization term
            delta = fmin(1.0 / delta_0, delta / delta_0);
            mu = mu * delta > mu_min ? mu * delta : 0.0;

            // Increment the number of backward passes
            num_backward_pass++;
        } while (num_backward_pass < max_num_backward_pass && !is_Q_uu_reg_positive_definite);

        // TEMP
        //std::cout << num_backward_pass << std::endl;

        // Increase the regularization term
        delta = fmax(delta_0, delta * delta_0);
        mu = fmax(mu_min, mu * delta);

        // Initialize counter for the number of forward passes
        int num_forward_pass = 0;

        // Allocate floats for the cost of the new trajectory and the expected cost reduction
        double J_new;
        double delta_J;

        // Boolean for whether or not to restart the forward pass
        bool is_restart_forward_pass;

        do {
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
                x_new[k + 1] = this->m_dynamics_ptr->f(x_new[k], u_new[k], t[k], dt);
            }

            // Compute the expected cost reduction of the control update
            delta_J = 0.0;
            for (int k = 0; k < this->m_num_discretization - 1; k++) {
                delta_J = delta_J + arma::as_scalar(this->m_learning_rate * gain_ff[k].t() * Q_u[k])
                        + arma::as_scalar(0.5 * pow(this->m_learning_rate, 2) * gain_ff[k].t() * Q_uu[k] * gain_ff[k]);
            }

            // Compute total cost of new trajectory
            J_new = this->m_cost_ptr->phi(x_new.back(), x_star);
            for (int k = 0; k < this->m_num_discretization - 1; k++) {
                J_new = J_new + this->m_cost_ptr->L(x_new[k], u_new[k], dt);
            }

            // Compute the cost reduction ratio
            double cost_reduction_ratio = (J_new - J) / delta_J;

            // Evalueate the forward pass
            if ((J_new - J) / delta_J < 0.1
                    || std::isinf(cost_reduction_ratio) || std::isnan(cost_reduction_ratio)) {
                // Decrease the learning rate
                // TEMP
                this->m_learning_rate = this->m_learning_rate * 0.5;

                // Restart the forward pass
                is_restart_forward_pass = true;
            } else {
                // Do not restart the forward pass
                is_restart_forward_pass = false;
            }

            // Increment number of forward passes
            num_forward_pass++;

            // TEMP
            //std::cout << this->m_learning_rate << "," << J_new << "," << J << "," << delta_J << "," << (J_new - J) / delta_J << std::endl;


        } while (num_forward_pass < max_num_forward_pass && is_restart_forward_pass);

        this->m_learning_rate = 1.0;

        // TEMP
        //std::cout << num_forward_pass << std::endl;

        // TEMP
        //std::cout << i << "," << J_new << std::endl;

        // Copy the updated control into the current nominal control
        u = u_new;

        // Copy the new trajectory into the old trajectory for the purpose of comparing trajectories between iterations
        x = x_new;
    }

    this->m_u_initial = u;

    // TEMP
    /**
    for (int k = 0; k < this->m_num_discretization; k++) {
        std::cout << x[k](0) << "," << x[k](1) << "," << x[k](2) << ",";
    }
    std::cout << "end" << std::endl;
     */

    // Return the first control input in the optimal control sequence
    return u[0];
}