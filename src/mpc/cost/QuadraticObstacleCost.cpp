/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "QuadraticObstacleCost.h"

double QuadraticObstacleCost::phi(arma::vec x_f, arma::vec x_star) {
    // Compute and return the quadratic terminal cost function
    return arma::as_scalar( 0.5 * (x_f - x_star).t() * this->m_Q_f * (x_f - x_star) );
}

arma::vec QuadraticObstacleCost::phi_x(arma::vec x_f, arma::vec x_star) {
    // Compute and return the gradient of the quadratic cost function
    return this->m_Q_f * (x_f - x_star);
}

arma::mat QuadraticObstacleCost::phi_xx(arma::vec x_f, arma::vec x_star) {
    // Compute and return the hessian of the quadratic cost function
    return this->m_Q_f;
}

double QuadraticObstacleCost::L(arma::vec x, arma::vec u, double dt) {
    // Declare transition cost and initialize with control cost
    double L = arma::as_scalar( 0.5 * u.t() * this->m_R * u * dt );

    // Compute and add the obstacle cost for each obstacle to the transition cost
    for (arma::vec obstacle : this->m_obstacles) {
        L = L + arma::as_scalar( exp(-0.5 * (x - obstacle).t() * this->m_sigma * (x - obstacle)) );
    }

    // Return the transition cost
    return L;
}

arma::vec QuadraticObstacleCost::L_x(arma::vec x, arma::vec u, double dt) {
    // Declare derivative of transition cost and initialize with zero
    arma::vec L_x = arma::zeros<arma::vec>(x.n_elem);

    // Compute derivative of running cost for each obstacle and add to the overall derivative
    for (arma::vec obstacle : this->m_obstacles) {
        L_x = L_x + arma::as_scalar( -exp(-0.5 * (x - obstacle).t() * this->m_sigma * (x - obstacle)) )
                        * this->m_sigma * (x - obstacle);
    }

    // Return the derivative of the transition cost
    return L_x;
}

arma::vec QuadraticObstacleCost::L_u(arma::vec x, arma::vec u, double dt) {
    // Compute and return the gradient of the transition cost with respect to the control
    return this->m_R * u * dt;
}

arma::mat QuadraticObstacleCost::L_xx(arma::vec x, arma::vec u, double dt) {
    // Declare second derivative of transition cost and initialize with zero
    arma::mat L_xx = arma::zeros<arma::mat>(x.n_elem, x.n_elem);

    // Compute second derivative of running cost for each obstacle and add to the overall derivative
    for (arma::vec obstacle : this->m_obstacles) {
        L_xx = L_xx + arma::as_scalar( exp(-0.5 * (x - obstacle).t() * this->m_sigma * (x - obstacle)) )
                         * (this->m_sigma * (x - obstacle) * (x - obstacle).t() * this->m_sigma - this->m_sigma);
    }

    // Return the derivative of the transition cost
    return L_xx;
}

arma::mat QuadraticObstacleCost::L_uu(arma::vec x, arma::vec u, double dt) {
    // Compute and return the Hessian of the transition cost with respect to the control
    return this->m_R * dt;
}

arma::mat QuadraticObstacleCost::L_xu(arma::vec x, arma::vec u, double dt) {
    // Compute and return the double double gradient of the transition cost with respect to the state and then control
    return arma::zeros<arma::mat>(x.n_elem, u.n_elem);
}

arma::mat QuadraticObstacleCost::L_ux(arma::vec x, arma::vec u, double dt) {
    // Compute and return the double double gradient of the transition cost with respect to the control and then state
    return arma::zeros<arma::mat>(u.n_elem, x.n_elem);
}