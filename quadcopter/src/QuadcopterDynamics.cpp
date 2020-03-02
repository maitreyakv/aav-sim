/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "QuadcopterDynamics.h"

arma::vec QuadcopterDynamics::F(arma::vec x, arma::vec u, double t) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = cos(x[6]);
    double x3 = sin(x[7]);
    double x4 = cos(x[8]);
    double x5 = x3*x4;
    double x6 = 1/this->m_m;
    double x7 = this->m_g*this->m_m;
    double x8 = this->m_k*(1000*u[0] + 1000*u[1] + 1000*u[2] + 1000*u[3] + x7/this->m_k);
    double x9 = x6*x8;
    double x10 = cos(x[7]);
    double x11 = x[11]*x10;
    double x12 = 1000*this->m_k*this->m_l;
    double x13 = x[10]*x2;
    double x14 = x0*x11 + x13;
    double x15 = x[10]*x0;
    double x16 = x11*x2 - x15;
    double x17 = x14*x16;
    double x18 = x[9]*x10;
    double x19 = -u[3];
    double x20 = x[9] - x[11]*x3;
    double x21 = x14*x20;
    double x22 = x[9]*x13 + x[11]*(x0*x18 + x13*x3) + (this->m_I_xx*x21 - this->m_I_yy*x21 + 1000*this->m_b*(u[0] - u[1] + u[2] + x19))/this->m_I_zz;
    double x23 = x2*x22;
    double x24 = 1/(pow(x0, 2) + pow(x2, 2));
    double x25 = x24/x10;
    double x26 = x25*x3;
    double x27 = x16*x20;
    double x28 = x[9]*x15 - x[11]*(-x15*x3 + x18*x2) + (-this->m_I_xx*x27 + this->m_I_zz*x27 + x12*(u[1] + x19))/this->m_I_yy;
    double x29 = x0*x28;

    // Compute linear acceleration components
    double X_ddot = x9*(x0*x1 + x2*x5);
    double Y_ddot = x9*(x0*x5 - x1*x2);
    double Z_ddot = -x6*(-x10*x4*x8 + x7);

    // Compute orientation angle "acceleration" components
    double phi_ddot = x[10]*x11 + x23*x26 + x26*x29 + (this->m_I_yy*x17 - this->m_I_zz*x17 + x12*(u[0] - u[2]))/this->m_I_xx;
    double theta_ddot = x24*(-x0*x22 + x2*x28);
    double psi_ddot = x25*(x23 + x29);

    // Allocate vector for the state derivative
    arma::vec dxdt = {x[3], x[4], x[5], X_ddot, Y_ddot, Z_ddot, x[9], x[10], x[11], phi_ddot, theta_ddot, psi_ddot};

    // Return the state derivative
    return dxdt;
}

arma::mat QuadcopterDynamics::Phi(arma::vec x, arma::vec u, double dt) {
    // Common subexpressions
    double x0 = sin(x[8]);
    double x1 = cos(x[6]);
    double x2 = x0*x1;
    double x3 = sin(x[7]);
    double x4 = sin(x[6]);
    double x5 = cos(x[8]);
    double x6 = x4*x5;
    double x7 = this->m_k*(this->m_g*this->m_m/this->m_k + 1000*u[0] + 1000*u[1] + 1000*u[2] + 1000*u[3])/this->m_m;
    double x8 = x1*x5;
    double x9 = cos(x[7]);
    double x10 = x7*x9;
    double x11 = x0*x4;
    double x12 = 1/this->m_I_yy;
    double x13 = x[11]*x3;
    double x14 = x[9] - x13;
    double x15 = x[10]*x1;
    double x16 = x[11]*x9;
    double x17 = x15 + x16*x4;
    double x18 = x14*x17;
    double x19 = x[9]*x9;
    double x20 = x15*x3 + x19*x4;
    double x21 = x[9]*x15 + x[11]*x20;
    double x22 = -x12*x18*(-this->m_I_xx + this->m_I_zz) + x21;
    double x23 = x22*x4;
    double x24 = pow(x4, 2);
    double x25 = pow(x1, 2);
    double x26 = 1/(x24 + x25);
    double x27 = 1/x9;
    double x28 = x26*x27;
    double x29 = x28*x3;
    double x30 = 1/this->m_I_xx;
    double x31 = pow(x17, 2);
    double x32 = x[10]*x4;
    double x33 = x1*x16;
    double x34 = -x32 + x33;
    double x35 = pow(x34, 2);
    double x36 = x14*x34;
    double x37 = 1/this->m_I_zz;
    double x38 = x37*(this->m_I_xx - this->m_I_yy);
    double x39 = -x36*x38;
    double x40 = x1*x19;
    double x41 = x3*x32;
    double x42 = x40 - x41;
    double x43 = x[9]*x32;
    double x44 = -x[11]*x42 + x43;
    double x45 = -u[3];
    double x46 = x12*(-this->m_I_xx*x36 + this->m_I_zz*x36 + 1000*this->m_k*this->m_l*(u[1] + x45)) + x44;
    double x47 = x1*x46;
    double x48 = x21 + x37*(this->m_I_xx*x18 - this->m_I_yy*x18 + 1000*this->m_b*(u[0] - u[1] + u[2] + x45));
    double x49 = x4*x48;
    double x50 = x1*x17;
    double x51 = x34*x4;
    double x52 = x30*(this->m_I_yy*x50 + this->m_I_yy*x51 - this->m_I_zz*x50 - this->m_I_zz*x51);
    double x53 = x1*x48;
    double x54 = x26*x53;
    double x55 = x4*x46;
    double x56 = x26*x55;
    double x57 = x14*x3;
    double x58 = x4*x57;
    double x59 = this->m_I_xx*x58;
    double x60 = this->m_I_yy*x58;
    double x61 = x17*x9;
    double x62 = this->m_I_xx*x61;
    double x63 = this->m_I_yy*x61;
    double x64 = x[9]*x3;
    double x65 = -x15*x9 + x4*x64;
    double x66 = x1*(x37*(x59 - x60 + x62 - x63) + x65);
    double x67 = x13*x28;
    double x68 = x1*x57;
    double x69 = x34*x9;
    double x70 = x1*x64 + x12*(this->m_I_xx*x68 + this->m_I_xx*x69 - this->m_I_zz*x68 - this->m_I_zz*x69) + x32*x9;
    double x71 = pow(x3, 2)/pow(x9, 2);
    double x72 = x17*x38 + x17;
    double x73 = x12*(this->m_I_xx - this->m_I_zz);
    double x74 = x34*x73;
    double x75 = x28*(x1*x72 - x4*(x34 + x74));
    double x76 = x14*x73;
    double x77 = x[9] + x13;
    double x78 = x24*(x76 + x77);
    double x79 = x14*x38;
    double x80 = x25*(x77 + x79);
    double x81 = x17*x4;
    double x82 = x1*x34;
    double x83 = x14*x9;
    double x84 = x4*x83;
    double x85 = x17*x3;
    double x86 = x20 + x37*(this->m_I_xx*x84 - this->m_I_xx*x85 - this->m_I_yy*x84 + this->m_I_yy*x85);
    double x87 = x1*x86;
    double x88 = x1*x83;
    double x89 = x3*x34;
    double x90 = x12*(this->m_I_xx*x88 - this->m_I_xx*x89 - this->m_I_zz*x88 + this->m_I_zz*x89);
    double x91 = x4*(x42 + x90);
    double x92 = -x40 + x41;
    double x93 = x[11]*x92 + x39 + x43;
    double x94 = x[11]*x70;
    double x95 = x27*x3;

    // Allocate matrix for Jacobian with respect to the state
    arma::mat F_x = arma::zeros<arma::mat>(x.n_elem, x.n_elem);

    // Compute components of Jacobian
    F_x(0,3) = 1;
    F_x(1,4) = 1;
    F_x(2,5) = 1;
    F_x(3,6) = -x7*(-x2 + x3*x6);
    F_x(3,7) = x10*x8;
    F_x(3,8) = x7*(-x2*x3 + x6);
    F_x(4,6) = x7*(x11 + x3*x8);
    F_x(4,7) = x10*x6;
    F_x(4,8) = -x7*(x11*x3 + x8);
    F_x(5,7) = -x3*x5*x7;
    F_x(5,8) = -x0*x10;
    F_x(6,9) = 1;
    F_x(7,10) = 1;
    F_x(8,11) = 1;
    F_x(9,6) = -x1*x29*(x39 + x44) + x23*x29 + x29*x47 - x29*x49 - x30*(this->m_I_yy*x31 - this->m_I_yy*x35 - this->m_I_zz*x31 + this->m_I_zz*x35);
    F_x(9,7) = -x[10]*x13 - x13*x52 + x4*x67*x70 + x54*x71 + x54 + x56*x71 + x56 - x66*x67;
    F_x(9,9) = x3*x75;
    F_x(9,10) = x16 + x29*x78 + x29*x80 - x30*(this->m_I_yy*x81 - this->m_I_yy*x82 - this->m_I_zz*x81 + this->m_I_zz*x82);
    F_x(9,11) = x[10]*x9 + x29*x87 - x29*x91 + x52*x9;
    F_x(10,6) = x26*(x1*x22 + x4*x93 - x53 - x55);
    F_x(10,7) = x26*(x[11]*x4*(-x37*(-x59 + x60 - x62 + x63) + x65) + x1*x94);
    F_x(10,9) = x26*(x1*(x32 - x33 - x74) - x4*x72);
    F_x(10,10) = x1*x26*x4*(x76 - x79);
    F_x(10,11) = x26*(x1*(-x90 + x92) - x4*x86);
    F_x(11,6) = x28*(-x1*x93 + x23 + x47 - x49);
    F_x(11,7) = x28*(-x[11]*x66 + x4*x94 + x53*x95 + x55*x95);
    F_x(11,9) = x75;
    F_x(11,10) = x28*(x78 + x80);
    F_x(11,11) = x28*(x87 - x91);

    // Transform the Jacobian into Phi and return
    return arma::eye<arma::mat>(x.n_elem, x.n_elem) + F_x * dt;
}

arma::mat QuadcopterDynamics::Beta(arma::vec x, arma::vec u, double dt) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = cos(x[6]);
    double x3 = cos(x[8]);
    double x4 = sin(x[7]);
    double x5 = x3*x4;
    double x6 = 1000*this->m_k/this->m_m;
    double x7 = x6*(x0*x1 + x2*x5);
    double x8 = x6*(x0*x5 - x1*x2);
    double x9 = cos(x[7]);
    double x10 = x3*x6*x9;
    double x11 = this->m_k*this->m_l;
    double x12 = x11/this->m_I_xx;
    double x13 = 1/(pow(x0, 2) + pow(x2, 2));
    double x14 = 1/x9;
    double x15 = this->m_b/this->m_I_zz;
    double x16 = x15*x2;
    double x17 = x13*x14*x16*x4;
    double x18 = x11/this->m_I_yy;
    double x19 = x0*x18;
    double x20 = 1000*x13;
    double x21 = x14*x20;
    double x22 = x21*(-x16 + x19);
    double x23 = x21*(x16 + x19);
    double x24 = x0*x15;
    double x25 = -x20*x24;
    double x26 = x18*x2;
    double x27 = x16*x21;

    // Allocate matrix for Jacobian with respect to the control
    arma::mat F_u = arma::zeros<arma::mat>(x.n_elem, u.n_elem);

    // Compute components of Jacobian
    F_u(3,0) = x7;
    F_u(3,1) = x7;
    F_u(3,2) = x7;
    F_u(3,3) = x7;
    F_u(4,0) = x8;
    F_u(4,1) = x8;
    F_u(4,2) = x8;
    F_u(4,3) = x8;
    F_u(5,0) = x10;
    F_u(5,1) = x10;
    F_u(5,2) = x10;
    F_u(5,3) = x10;
    F_u(9,0) = 1000*x12 + 1000*x17;
    F_u(9,1) = x22*x4;
    F_u(9,2) = -1000*x12 + 1000*x17;
    F_u(9,3) = -x23*x4;
    F_u(10,0) = x25;
    F_u(10,1) = x20*(x24 + x26);
    F_u(10,2) = x25;
    F_u(10,3) = x20*(x24 - x26);
    F_u(11,0) = x27;
    F_u(11,1) = x22;
    F_u(11,2) = x27;
    F_u(11,3) = -x23;

    // Transform the Jacobian into Beta and return
    return F_u * dt;
}
