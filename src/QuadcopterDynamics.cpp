/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "QuadcopterDynamics.h"

arma::vec QuadcopterDynamics::F(arma::vec x, arma::vec u, double t) {
    // Compute trig functions of orientation angles
    double s_phi   = sin(x[6]);
    double s_theta = sin(x[7]);
    double s_psi   = sin(x[8]);
    double c_phi   = cos(x[6]);
    double c_theta = cos(x[7]);
    double c_psi   = cos(x[8]);

    // Compute square of trig functions
    double s_phi_sq = pow(s_phi, 2);
    double s_theta_sq = pow(s_theta, 2);
    double s_psi_sq = pow(s_psi, 2);
    double c_phi_sq = pow(c_phi, 2);
    double c_theta_sq = pow(c_theta, 2);
    double c_psi_sq = pow(c_psi, 2);

    // Compute linear acceleration components
    double X_ddot = this->m_k*(u[0] + u[1] + u[2] + u[3])*s_psi*s_theta/this->m_m;
    double Y_ddot = -this->m_k*(u[0] + u[1] + u[2] + u[3])*s_theta*c_psi/this->m_m;
    double Z_ddot = (-this->m_g*this->m_m + this->m_k*(u[0] + u[1] + u[2] + u[3])*c_theta)/this->m_m;

    // Compute orientation angle "acceleration" components
    double phi_ddot = x[11]*x[10]*c_theta + (x[9]*x[10]*s_phi - x[11]*(x[9]*c_phi*c_theta - x[10]*s_phi*s_theta)
                    + (-this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi) + this->m_I_zz*(x[9]
                    - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                    + this->m_kl*(u[1] - u[3]))/this->m_I_yy)*s_phi*s_theta/(s_phi_sq*c_theta + c_phi_sq*c_theta)
                    + (x[9]*x[10]*c_phi - x[11]*(-x[9]*s_phi*c_theta - x[10]*s_theta*c_phi)
                    + (this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi)
                    - this->m_I_yy*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi) + this->m_b*(u[0] - u[1]
                    + u[2] - u[3]))/this->m_I_zz)*s_theta*c_phi/(s_phi_sq*c_theta + c_phi_sq*c_theta)
                    + (this->m_I_yy*(x[11]*s_phi*c_theta + x[10]*c_phi)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                    - this->m_I_zz*(x[11]*s_phi*c_theta + x[10]*c_phi)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                    + this->m_kl*(u[0] - u[2]))/this->m_I_xx;
    double theta_ddot = (x[9]*x[10]*s_phi - x[11]*(x[9]*c_phi*c_theta - x[10]*s_phi*s_theta)
                      + (-this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                      + this->m_I_zz*(x[9] - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                      + this->m_kl*(u[1] - u[3]))/this->m_I_yy)*c_phi*c_theta/(s_phi_sq*c_theta + c_phi_sq*c_theta)
                      - (x[9]*x[10]*c_phi - x[11]*(-x[9]*s_phi*c_theta - x[10]*s_theta*c_phi)
                      + (this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi)
                      - this->m_I_yy*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi) + this->m_b*(u[0]
                      - u[1] + u[2] - u[3]))/this->m_I_zz)*s_phi*c_theta/(s_phi_sq*c_theta + c_phi_sq*c_theta);
    double psi_ddot = (x[9]*x[10]*s_phi - x[11]*(x[9]*c_phi*c_theta - x[10]*s_phi*s_theta)
                    + (-this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                    + this->m_I_zz*(x[9] - x[11]*s_theta)*(x[11]*c_phi*c_theta - x[10]*s_phi)
                    + this->m_kl*(u[1] - u[3]))/this->m_I_yy)*s_phi/(s_phi_sq*c_theta + c_phi_sq*c_theta)
                    + (x[9]*x[10]*c_phi - x[11]*(-x[9]*s_phi*c_theta - x[10]*s_theta*c_phi)
                    + (this->m_I_xx*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi)
                    - this->m_I_yy*(x[9] - x[11]*s_theta)*(x[11]*s_phi*c_theta + x[10]*c_phi)
                    + this->m_b*(u[0] - u[1] + u[2] - u[3]))/this->m_I_zz)*c_phi/(s_phi_sq*c_theta + c_phi_sq*c_theta);

    // Allocate vector for the state derivative
    arma::vec dxdt = {x[3], x[4], x[5], X_ddot, Y_ddot, Z_ddot, x[9], x[10], x[11], phi_ddot, theta_ddot, psi_ddot};

    // Return the state derivative
    return dxdt;
}

arma::mat QuadcopterDynamics::Phi(arma::vec x, arma::vec u, double dt) {
    // TODO: IMPLEMENT ME
}

arma::mat QuadcopterDynamics::Beta(arma::vec x, arma::vec u, double dt) {
    // TODO: IMPLEMENT ME
}