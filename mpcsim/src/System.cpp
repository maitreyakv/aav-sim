/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "System.h"

void System::operator()(const arma::vec &x, arma::vec &dxdt, const double t) {
    // Use the Dynamics of the System to compute F using the current state and the current input
    dxdt = this->m_dynamics_ptr->F(x, this->m_u, t);
}

void System::updateControl(arma::vec x, arma::vec x_star, double t, double horizon) {
    // Use the controller to compute the best control input for the current state and assign it to the member variable
    this->m_u = this->m_controller_ptr->computeOptimalControl(x, t, x_star, t + horizon);

#ifdef DEBUG
    for (long long unsigned int i = 0; i < this->m_u.n_elem; i++) {
        if (std::isnan(this->m_u(i))) {
            std::cout << "error: computed optimal control has NaN\n"; exit(-1);
        }
    }
#endif
}