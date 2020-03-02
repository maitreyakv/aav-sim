/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "QuadraticCost.h"

double QuadraticCost::phi(arma::vec x_f, arma::vec x_star) {
    // Compute and return the quadratic terminal cost function
    return arma::as_scalar( 0.5 * (x_f - x_star).t() * this->m_Q_f * (x_f - x_star) );
}

arma::vec QuadraticCost::phi_x(arma::vec x_f, arma::vec x_star) {
    // Compute and return the gradient of the quadratic cost function
    return this->m_Q_f * (x_f - x_star);
}

arma::mat QuadraticCost::phi_xx(arma::vec x_f, arma::vec x_star) {
    // Compute and return the hessian of the quadratic cost function
    return this->m_Q_f;
}

double QuadraticCost::L(arma::vec x, arma::vec u, double dt) {
    // Compute and return the transition cost
    return arma::as_scalar( 0.5 * u.t() * this->m_R * u * dt );
}

arma::vec QuadraticCost::L_x(arma::vec x, arma::vec u, double dt) {
    // Compute and return the gradient of the transition cost with respect to the state
    return arma::zeros<arma::vec>(x.n_elem);
}

arma::vec QuadraticCost::L_u(arma::vec x, arma::vec u, double dt) {
    // Compute and return the gradient of the transition cost with respect to the control
    return this->m_R * u * dt;
}

arma::mat QuadraticCost::L_xx(arma::vec x, arma::vec u, double dt) {
    // Compute and return the Hessian of the transition cost with respect to the state
    return arma::zeros<arma::mat>(x.n_elem, x.n_elem);
}

arma::mat QuadraticCost::L_uu(arma::vec x, arma::vec u, double dt) {
    // Compute and return the Hessian of the transition cost with respect to the control
    return this->m_R * dt;
}

arma::mat QuadraticCost::L_xu(arma::vec x, arma::vec u, double dt) {
    // Compute and return the double double gradient of the transition cost with respect to the state and then control
    return arma::zeros<arma::mat>(x.n_elem, u.n_elem);
}

arma::mat QuadraticCost::L_ux(arma::vec x, arma::vec u, double dt) {
    // Compute and return the double double gradient of the transition cost with respect to the control and then state
    return arma::zeros<arma::mat>(u.n_elem, x.n_elem);
}