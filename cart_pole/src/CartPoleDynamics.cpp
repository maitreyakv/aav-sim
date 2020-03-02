/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "CartPoleDynamics.h"

arma::vec CartPoleDynamics::F(arma::vec x, arma::vec u, double t) {
    // Compute sine and cosine of pendulum angle theta, for use in multiple places
    double sin_theta = sin( x(2) );
    double cos_theta = cos( x(2) );

    // Compute intermediate terms in the dynamics
    double lambda = this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2);
    double psi = u(0) + this->m_mass_pole * sin_theta * (this->m_length_pole * pow(x(3), 2) + 9.81 * cos_theta);

    double eta = -u(0) * cos_theta - this->m_mass_pole * this->m_length_pole * pow(x(3), 2) * cos_theta * sin_theta
                 - (this->m_mass_cart + this->m_mass_pole) * 9.81 * sin_theta;

    // Compute the cart acceleration
    double z_ddot = psi / lambda;

    // Compute the pole angular acceleration
    double theta_ddot = eta / (this->m_length_pole * lambda);

    // Assemble the state derivative
    arma::vec dxdt = { x(1), z_ddot, x(3), theta_ddot };

    // Return the state derivative
    return dxdt;
}

arma::mat CartPoleDynamics::Phi(arma::vec x, arma::vec u, double dt) {
    // Compute sine and cosine of pendulum angle theta, for use in multiple places
    double sin_theta = sin( x(2) );
    double cos_theta = cos( x(2) );

    // Allocate a vector for the gradient of F with respect to the state
    arma::mat F_x = arma::zeros<arma::mat>(x.n_elem, x.n_elem);

    // Compute terms in the matrix Jacobian matrix
    F_x(0,1) = 1.0;
    F_x(1,2) = (this->m_mass_pole * cos_theta * (this->m_length_pole * pow(x(3), 2) + 9.81 * cos_theta) - 9.81
               * this->m_mass_pole * pow(sin_theta, 2)) / (this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2))
               - (2.0 * this->m_mass_pole * cos_theta * sin_theta * (u(0) + this->m_mass_pole * sin_theta
               * (this->m_length_pole * pow(x(3), 2) + 9.81 * cos_theta)))
               / pow(this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2), 2);
    F_x(1,3) = (2.0 * this->m_length_pole * this->m_mass_pole * x(3) * sin_theta)
               / (this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2));
    F_x(2,3) = 1.0;
    F_x(3,2) = (u(0) * sin_theta - 9.81 * cos_theta * this->m_mass_total - this->m_length_pole
               * this->m_mass_pole * pow(x(3), 2) * pow(cos_theta, 2) + this->m_length_pole * this->m_mass_pole
               * pow(x(3), 2) * pow(sin_theta, 2)) / (this->m_length_pole * (this->m_mass_cart + this->m_mass_pole
               * pow(sin_theta, 2))) + (2.0 * this->m_mass_pole * cos_theta * sin_theta * (this->m_length_pole
               * this->m_mass_pole * cos_theta * sin_theta * pow(x(3), 2) + u(0) * cos_theta + 9.81 * sin_theta
               * this->m_mass_total))
               / (this->m_length_pole * pow(this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2), 2));
    F_x(3,3) = -(2.0 * this->m_mass_pole * x(3) * cos_theta * sin_theta)
               / (this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2));

    // Transform the Jacobian into Phi and return
    return arma::eye<arma::mat>(x.n_elem, x.n_elem) + F_x * dt;
}

arma::mat CartPoleDynamics::Beta(arma::vec x, arma::vec u, double dt) {
    // Compute sine and cosine of pendulum angle theta, for use in multiple places
    double sin_theta = sin( x(2) );
    double cos_theta = cos( x(2) );

    // Allocate a vector for the gradient of F with respect to the control
    arma::mat F_u = arma::zeros<arma::vec>(x.n_elem, 1);

    // Compute and add the components of the gradient to the vector
    F_u(1) = 1.0 / ( this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2) );
    F_u(3) = -cos_theta / ( this->m_length_pole * ( this->m_mass_cart + this->m_mass_pole * pow(sin_theta, 2) ) );

    // Compute and return the scaled gradient, beta
    return F_u * dt;
}