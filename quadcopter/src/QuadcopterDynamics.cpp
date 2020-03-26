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
    double x8 = sqrt(x7/this->m_k);
    double x9 = pow(2*u[1] + x8, 2);
    double x10 = pow(2*u[3] + x8, 2);
    double x11 = pow(2*u[0] + x8, 2);
    double x12 = pow(2*u[2] + x8, 2);
    double x13 = x11 + x12;
    double x14 = this->m_k*(x10 + x13 + x9)/4;
    double x15 = x14*x6;
    double x16 = cos(x[7]);
    double x17 = x[11]*x16;
    double x18 = x[10]*x2;
    double x19 = x0*x17 + x18;
    double x20 = x[10]*x0;
    double x21 = 4*x17*x2 - 4*x20;
    double x22 = x19*x21;
    double x23 = this->m_k*this->m_l;
    double x24 = 4*x[9];
    double x25 = x[9]*x16;
    double x26 = 4*x[11];
    double x27 = x[9] - x[11]*x3;
    double x28 = x21*x27;
    double x29 = -x10;
    double x30 = x20*x24 - x26*(x2*x25 - x20*x3) + (-this->m_I_xx*x28 + this->m_I_zz*x28 + x23*(x29 + x9))/this->m_I_yy;
    double x31 = x0*x30;
    double x32 = 1/(4*(pow(x0, 2) + pow(x2, 2)));
    double x33 = x32/x16;
    double x34 = x3*x33;
    double x35 = 4*x19*x27;
    double x36 = x18*x24 + x26*(x0*x25 + x18*x3) + (this->m_I_xx*x35 - this->m_I_yy*x35 + this->m_b*(x13 + x29 - x9))/this->m_I_zz;
    double x37 = x2*x36;

    // Compute linear acceleration components
    double X_ddot = x15*(x0*x1 + x2*x5);
    double Y_ddot = x15*(x0*x5 - x1*x2);
    double Z_ddot = -x6*(-x14*x16*x4 + x7);

    // Compute orientation angle "acceleration" components
    double phi_ddot = x[10]*x17 + x31*x34 + x34*x37 + (this->m_I_yy*x22 - this->m_I_zz*x22 + x23*(x11 - x12))/(4*this->m_I_xx);
    double theta_ddot = x32*(-x0*x36 + x2*x30);
    double psi_ddot = x33*(x31 + x37);

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
    double x7 = sqrt(this->m_g*this->m_m/this->m_k);
    double x8 = pow(2*u[1] + x7, 2);
    double x9 = pow(2*u[3] + x7, 2);
    double x10 = pow(2*u[0] + x7, 2) + pow(2*u[2] + x7, 2);
    double x11 = this->m_k*(x10 + x8 + x9)/(4*this->m_m);
    double x12 = x1*x5;
    double x13 = cos(x[7]);
    double x14 = x11*x13;
    double x15 = x0*x4;
    double x16 = x[10]*x1;
    double x17 = x[9]*x16;
    double x18 = x[9]*x13;
    double x19 = x16*x3 + x18*x4;
    double x20 = x[11]*x19;
    double x21 = 1/this->m_I_yy;
    double x22 = x[11]*x13;
    double x23 = x16 + x22*x4;
    double x24 = x[11]*x3;
    double x25 = x[9] - x24;
    double x26 = x23*x25;
    double x27 = x17 + x20 - x21*x26*(-this->m_I_xx + this->m_I_zz);
    double x28 = x27*x4;
    double x29 = pow(x4, 2);
    double x30 = pow(x1, 2);
    double x31 = 1/(x29 + x30);
    double x32 = 1/x13;
    double x33 = x31*x32;
    double x34 = x3*x33;
    double x35 = 1/this->m_I_xx;
    double x36 = pow(x23, 2);
    double x37 = x[10]*x4;
    double x38 = x1*x22;
    double x39 = -x37 + x38;
    double x40 = pow(x39, 2);
    double x41 = x1*x18;
    double x42 = x3*x37;
    double x43 = x41 - x42;
    double x44 = x[11]*x43;
    double x45 = x[9]*x37;
    double x46 = 1/this->m_I_zz;
    double x47 = x46*(this->m_I_xx - this->m_I_yy);
    double x48 = x25*x47;
    double x49 = -x39*x48 + x45;
    double x50 = x1*x3;
    double x51 = this->m_I_xx*x25;
    double x52 = 4*x39;
    double x53 = this->m_I_zz*x25;
    double x54 = -x9;
    double x55 = x21*(this->m_k*this->m_l*(x54 + x8) - x51*x52 + x52*x53) - 4*x44 + 4*x45;
    double x56 = x1/4;
    double x57 = x55*x56;
    double x58 = 4*x26;
    double x59 = 4*x17 + 4*x20 + x46*(this->m_I_xx*x58 - this->m_I_yy*x58 + this->m_b*(x10 + x54 - x8));
    double x60 = x4/4;
    double x61 = x59*x60;
    double x62 = x1*x23;
    double x63 = x39*x4;
    double x64 = x35*(this->m_I_yy*x62 + this->m_I_yy*x63 - this->m_I_zz*x62 - this->m_I_zz*x63);
    double x65 = x3*x4;
    double x66 = x51*x65;
    double x67 = this->m_I_yy*x25;
    double x68 = x65*x67;
    double x69 = x13*x23;
    double x70 = this->m_I_xx*x69;
    double x71 = this->m_I_yy*x69;
    double x72 = x[9]*x3;
    double x73 = -x13*x16 + x4*x72;
    double x74 = x1*(x46*(x66 - x68 + x70 - x71) + x73);
    double x75 = x24*x33;
    double x76 = x13*x39;
    double x77 = x1*x72 + x13*x37 + x21*(this->m_I_xx*x76 - this->m_I_zz*x76 + x50*x51 - x50*x53);
    double x78 = x55*x60;
    double x79 = x31*x78;
    double x80 = pow(x3, 2)/pow(x13, 2);
    double x81 = x56*x59;
    double x82 = x31*x81;
    double x83 = x23*x47 + x23;
    double x84 = x21*(this->m_I_xx - this->m_I_zz);
    double x85 = x39*x84;
    double x86 = x33*(x1*x83 - x4*(x39 + x85));
    double x87 = x25*x84;
    double x88 = x[9] + x24;
    double x89 = x29*(x87 + x88);
    double x90 = x30*(x48 + x88);
    double x91 = x23*x4;
    double x92 = x1*x39;
    double x93 = x13*x4;
    double x94 = x23*x3;
    double x95 = x19 + x46*(-this->m_I_xx*x94 + this->m_I_yy*x94 + x51*x93 - x67*x93);
    double x96 = x1*x95;
    double x97 = x1*x13;
    double x98 = x3*x39;
    double x99 = x21*(-this->m_I_xx*x98 + this->m_I_zz*x98 + x51*x97 - x53*x97);
    double x100 = x4*(x43 + x99);
    double x101 = -x41 + x42;
    double x102 = x[11]*x101 + x49;
    double x103 = x[11]*x77;
    double x104 = x3*x32;

    // Allocate matrix for Jacobian with respect to the state
    arma::mat F_x = arma::zeros<arma::mat>(x.n_elem, x.n_elem);

    // Compute components of Jacobian
    F_x(0,3) = 1;
    F_x(1,4) = 1;
    F_x(2,5) = 1;
    F_x(3,6) = -x11*(-x2 + x3*x6);
    F_x(3,7) = x12*x14;
    F_x(3,8) = x11*(-x2*x3 + x6);
    F_x(4,6) = x11*(x12*x3 + x15);
    F_x(4,7) = x14*x6;
    F_x(4,8) = -x11*(x12 + x15*x3);
    F_x(5,7) = -x11*x3*x5;
    F_x(5,8) = -x0*x14;
    F_x(6,9) = 1;
    F_x(7,10) = 1;
    F_x(8,11) = 1;
    F_x(9,6) = x28*x34 - x33*x50*(-x44 + x49) + x34*x57 - x34*x61 - x35*(this->m_I_yy*x36 - this->m_I_yy*x40 - this->m_I_zz*x36 + this->m_I_zz*x40);
    F_x(9,7) = -x[10]*x24 - x24*x64 + x4*x75*x77 - x74*x75 + x79*x80 + x79 + x80*x82 + x82;
    F_x(9,9) = x3*x86;
    F_x(9,10) = x22 + x34*x89 + x34*x90 - x35*(this->m_I_yy*x91 - this->m_I_yy*x92 - this->m_I_zz*x91 + this->m_I_zz*x92);
    F_x(9,11) = x[10]*x13 - x100*x34 + x13*x64 + x34*x96;
    F_x(10,6) = x31*(x1*x27 + x102*x4 - x78 - x81);
    F_x(10,7) = x31*(x[11]*x4*(-x46*(-x66 + x68 - x70 + x71) + x73) + x1*x103);
    F_x(10,9) = x31*(x1*(x37 - x38 - x85) - x4*x83);
    F_x(10,10) = x1*x31*x4*(-x48 + x87);
    F_x(10,11) = x31*(x1*(x101 - x99) - x4*x95);
    F_x(11,6) = x33*(-x1*x102 + x28 + x57 - x61);
    F_x(11,7) = x33*(-x[11]*x74 + x103*x4 + x104*x78 + x104*x81);
    F_x(11,9) = x86;
    F_x(11,10) = x33*(x89 + x90);
    F_x(11,11) = x33*(-x100 + x96);

    // Transform the Jacobian into Phi and return
    return arma::eye<arma::mat>(x.n_elem, x.n_elem) + F_x * dt;
}

arma::mat QuadcopterDynamics::Beta(arma::vec x, arma::vec u, double dt) {
    // Common subexpressions
    double x0 = sqrt(this->m_g*this->m_m/this->m_k);
    double x1 = 2*u[0] + x0;
    double x2 = sin(x[6]);
    double x3 = sin(x[8]);
    double x4 = cos(x[6]);
    double x5 = cos(x[8]);
    double x6 = sin(x[7]);
    double x7 = x5*x6;
    double x8 = this->m_k/this->m_m;
    double x9 = x8*(x2*x3 + x4*x7);
    double x10 = 2*u[1] + x0;
    double x11 = 2*u[2] + x0;
    double x12 = 2*u[3] + x0;
    double x13 = x8*(x2*x7 - x3*x4);
    double x14 = cos(x[7]);
    double x15 = x14*x5*x8;
    double x16 = this->m_k*this->m_l;
    double x17 = x16/this->m_I_xx;
    double x18 = 1/(pow(x2, 2) + pow(x4, 2));
    double x19 = 1/x14;
    double x20 = this->m_b/this->m_I_zz;
    double x21 = x20*x4;
    double x22 = x18*x19*x21;
    double x23 = x22*x6;
    double x24 = x16/this->m_I_yy;
    double x25 = x2*x24;
    double x26 = x10*x18;
    double x27 = x19*x26*(-x21 + x25);
    double x28 = x12*x18;
    double x29 = x19*x28*(x21 + x25);
    double x30 = x2*x20;
    double x31 = x18*x30;
    double x32 = x24*x4;

    // Allocate matrix for Jacobian with respect to the control
    arma::mat F_u = arma::zeros<arma::mat>(x.n_elem, u.n_elem);

    // Compute components of Jacobian
    F_u(3,0) = x1*x9;
    F_u(3,1) = x10*x9;
    F_u(3,2) = x11*x9;
    F_u(3,3) = x12*x9;
    F_u(4,0) = x1*x13;
    F_u(4,1) = x10*x13;
    F_u(4,2) = x11*x13;
    F_u(4,3) = x12*x13;
    F_u(5,0) = x1*x15;
    F_u(5,1) = x10*x15;
    F_u(5,2) = x11*x15;
    F_u(5,3) = x12*x15;
    F_u(9,0) = x1*(x17 + x23);
    F_u(9,1) = x27*x6;
    F_u(9,2) = x11*(-x17 + x23);
    F_u(9,3) = -x29*x6;
    F_u(10,0) = -x1*x31;
    F_u(10,1) = x26*(x30 + x32);
    F_u(10,2) = -x11*x31;
    F_u(10,3) = x28*(x30 - x32);
    F_u(11,0) = x1*x22;
    F_u(11,1) = x27;
    F_u(11,2) = x11*x22;
    F_u(11,3) = -x29;

    // Transform the Jacobian into Beta and return
    return F_u * dt;
}
