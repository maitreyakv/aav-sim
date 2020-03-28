/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "QuadcopterDynamics.h"

// TODO: Add documentation

arma::vec QuadcopterDynamics::F(arma::vec x, arma::vec u, double t) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = cos(x[6]);
    double x3 = sin(x[7]);
    double x4 = cos(x[8]);
    double x5 = x3*x4;
    double x6 = 1.0/this->m_m;
    double x7 = this->m_g*this->m_m;
    double x8 = sqrt(x7/this->m_k);
    double x9 = pow(2*u[1] + x8, 2);
    double x10 = pow(2*u[3] + x8, 2);
    double x11 = pow(2*u[0] + x8, 2);
    double x12 = pow(2*u[2] + x8, 2);
    double x13 = x11 + x12;
    double x14 = (1.0/4.0)*this->m_k*(x10 + x13 + x9);
    double x15 = x14*x6;
    double x16 = cos(x[7]);
    double x17 = x16*x[11];
    double x18 = x2*x[10];
    double x19 = x0*x17 + x18;
    double x20 = x0*x[10];
    double x21 = -4*x17*x2 + 4*x20;
    double x22 = x19*x21;
    double x23 = this->m_k*this->m_l;
    double x24 = 4*x[9];
    double x25 = x16*x[9];
    double x26 = 4*x[11];
    double x27 = x3*x[11] - x[9];
    double x28 = x21*x27;
    double x29 = -x10;
    double x30 = x20*x24 + x26*(-x2*x25 + x20*x3) + (-this->m_I_xx*x28 + this->m_I_zz*x28 + x23*(x29 + x9))/this->m_I_yy;
    double x31 = x0*x30;
    double x32 = (1.0/4.0)/(pow(x0, 2) + pow(x2, 2));
    double x33 = x32/x16;
    double x34 = x3*x33;
    double x35 = 4*x19*x27;
    double x36 = x18*x24 + x26*(x0*x25 + x18*x3) + (-this->m_I_xx*x35 + this->m_I_yy*x35 + this->m_b*(x13 + x29 - x9))/this->m_I_zz;
    double x37 = x2*x36;

    // Allocate vector for the function F
    arma::vec F = arma::zeros<arma::vec>(x.n_elem);

    // Compute components of function f
    F[0]  = x[3];
    F[1]  = x[4];
    F[2]  = x[5];
    F[3]  = x15*(x0*x1 + x2*x5);
    F[4]  = x15*(x0*x5 - x1*x2);
    F[5]  = -x6*(-x14*x16*x4 + x7);
    F[6]  = x[9];
    F[7]  = x[10];
    F[8]  = x[11];
    F[9]  = x17*x[10] + x31*x34 + x34*x37 + (1.0/4.0)*(-this->m_I_yy*x22 + this->m_I_zz*x22 + x23*(x11 - x12))/this->m_I_xx;
    F[10] = x32*(-x0*x36 + x2*x30);
    F[11] = x33*(x31 + x37);

    // Return the state derivative F
    return F;
}

arma::vec QuadcopterDynamics::f(arma::vec x, arma::vec u, double t, double dt) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = cos(x[6]);
    double x3 = sin(x[7]);
    double x4 = cos(x[8]);
    double x5 = x3*x4;
    double x6 = this->m_g*this->m_m;
    double x7 = sqrt(x6/this->m_k);
    double x8 = pow(2*u[1] + x7, 2);
    double x9 = pow(2*u[3] + x7, 2);
    double x10 = pow(2*u[0] + x7, 2);
    double x11 = pow(2*u[2] + x7, 2);
    double x12 = x10 + x11;
    double x13 = this->m_k*(x12 + x8 + x9);
    double x14 = (1.0/4.0)*dt;
    double x15 = x14/this->m_m;
    double x16 = x13*x15;
    double x17 = cos(x[7]);
    double x18 = x17*x[11];
    double x19 = x2*x[10];
    double x20 = x0*x18 + x19;
    double x21 = x0*x[10];
    double x22 = -4*x18*x2 + 4*x21;
    double x23 = x20*x22;
    double x24 = this->m_k*this->m_l;
    double x25 = 4*x[9];
    double x26 = x17*x[9];
    double x27 = 4*x[11];
    double x28 = x3*x[11] - x[9];
    double x29 = x22*x28;
    double x30 = -x9;
    double x31 = x21*x25 + x27*(-x2*x26 + x21*x3) + (-this->m_I_xx*x29 + this->m_I_zz*x29 + x24*(x30 + x8))/this->m_I_yy;
    double x32 = x0*x31;
    double x33 = 1.0/(pow(x0, 2) + pow(x2, 2));
    double x34 = 1.0/x17;
    double x35 = x3*x33*x34;
    double x36 = 4*x20*x28;
    double x37 = x19*x25 + x27*(x0*x26 + x19*x3) + (-this->m_I_xx*x36 + this->m_I_yy*x36 + this->m_b*(x12 + x30 - x8))/this->m_I_zz;
    double x38 = x2*x37;
    double x39 = x14*x33;

    // Allocate vector for the function f
    arma::vec f = arma::zeros<arma::vec>(x.n_elem);

    // Compute components of function f
    f[0]  = dt*x[3] + x[0];
    f[1]  = dt*x[4] + x[1];
    f[2]  = dt*x[5] + x[2];
    f[3]  = x16*(x0*x1 + x2*x5) + x[3];
    f[4]  = x16*(x0*x5 - x1*x2) + x[4];
    f[5]  = -x15*(-x13*x17*x4 + 4*x6) + x[5];
    f[6]  = dt*x[9] + x[6];
    f[7]  = dt*x[10] + x[7];
    f[8]  = dt*x[11] + x[8];
    f[9]  = x14*(4*x18*x[10] + x32*x35 + x35*x38 + (-this->m_I_yy*x23 + this->m_I_zz*x23 + x24*(x10 - x11))/this->m_I_xx) + x[9];
    f[10] = x39*(-x0*x37 + x2*x31) + x[10];
    f[11] = x34*x39*(x32 + x38) + x[11];

    // Return the next state determined by the equations of motion
    return f;
}

arma::mat QuadcopterDynamics::f_x(arma::vec x, arma::vec u, double t, double dt) {
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
    double x11 = (1.0/4.0)*dt*this->m_k*(x10 + x8 + x9)/this->m_m;
    double x12 = x1*x5;
    double x13 = cos(x[7]);
    double x14 = x11*x13;
    double x15 = x0*x4;
    double x16 = 1.0/this->m_I_xx;
    double x17 = x1*x[10];
    double x18 = x13*x[11];
    double x19 = x17 + x18*x4;
    double x20 = pow(x19, 2);
    double x21 = x4*x[10];
    double x22 = x1*x18;
    double x23 = x21 - x22;
    double x24 = pow(x23, 2);
    double x25 = x3*x[11];
    double x26 = x25 - x[9];
    double x27 = 1.0/this->m_I_yy;
    double x28 = x27*(this->m_I_xx - this->m_I_zz);
    double x29 = x26*x28;
    double x30 = x13*x[9];
    double x31 = x17*x3 + x30*x4;
    double x32 = x31*x[11];
    double x33 = x17*x[9];
    double x34 = x32 + x33;
    double x35 = x4*(-x19*x29 + x34);
    double x36 = pow(x4, 2);
    double x37 = pow(x1, 2);
    double x38 = 1.0/(x36 + x37);
    double x39 = 1.0/x13;
    double x40 = x3*x39;
    double x41 = x38*x40;
    double x42 = 1.0/this->m_I_zz;
    double x43 = x42*(this->m_I_xx - this->m_I_yy);
    double x44 = x26*x43;
    double x45 = -x1*x30 + x21*x3;
    double x46 = x45*x[11];
    double x47 = x21*x[9];
    double x48 = x46 + x47;
    double x49 = x1*(-x23*x44 + x48);
    double x50 = 4*x26;
    double x51 = x23*x50;
    double x52 = -x9;
    double x53 = x27*(-this->m_I_xx*x51 + this->m_I_zz*x51 + this->m_k*this->m_l*(x52 + x8)) + 4*x46 + 4*x47;
    double x54 = (1.0/4.0)*x1;
    double x55 = x53*x54;
    double x56 = x19*x50;
    double x57 = 4*x32 + 4*x33 + x42*(-this->m_I_xx*x56 + this->m_I_yy*x56 + this->m_b*(x10 + x52 - x8));
    double x58 = (1.0/4.0)*x4;
    double x59 = x57*x58;
    double x60 = x1*x19;
    double x61 = this->m_I_yy*x60;
    double x62 = x23*x4;
    double x63 = this->m_I_zz*x62;
    double x64 = this->m_I_zz*x60;
    double x65 = this->m_I_yy*x62;
    double x66 = x3*x[9];
    double x67 = this->m_I_xx*x3;
    double x68 = -x25 + x[9];
    double x69 = x4*x68;
    double x70 = x13*x19;
    double x71 = -x13*x17 + x4*x66 + x42*(this->m_I_xx*x70 - this->m_I_yy*x3*x69 - this->m_I_yy*x70 + x67*x69);
    double x72 = x25*x38*x39;
    double x73 = x1*x26;
    double x74 = this->m_I_zz*x73;
    double x75 = x13*x23;
    double x76 = x1*x66 + x13*x21 - x27*(this->m_I_xx*x75 - this->m_I_zz*x75 - x3*x74 + x67*x73);
    double x77 = x53*x58;
    double x78 = x38*x77;
    double x79 = pow(x3, 2)/pow(x13, 2);
    double x80 = x54*x57;
    double x81 = x38*x80;
    double x82 = x19*x43 + x19;
    double x83 = -x21 + x22;
    double x84 = x28*x83 + x83;
    double x85 = dt*x38;
    double x86 = x39*x85;
    double x87 = x86*(x1*x82 - x4*x84);
    double x88 = -x29;
    double x89 = x25 + x[9];
    double x90 = x36*(x88 + x89);
    double x91 = x37*(-x44 + x89);
    double x92 = x19*x4;
    double x93 = x1*x23;
    double x94 = this->m_I_xx*x13;
    double x95 = x26*x4;
    double x96 = x19*x3;
    double x97 = x31 - x42*(this->m_I_xx*x96 - this->m_I_yy*x13*x95 - this->m_I_yy*x96 + x94*x95);
    double x98 = x1*x97;
    double x99 = x23*x3;
    double x100 = -x27*(this->m_I_xx*x99 - this->m_I_zz*x99 + x13*x74 - x73*x94) + x45;
    double x101 = x100*x4;
    double x102 = x4*x[11];
    double x103 = x1*x[11];

    // Allocate matrix for Jacobian with respect to the state
    arma::mat f_x = arma::zeros<arma::mat>(x.n_elem, x.n_elem);

    // Compute components of Jacobian
    f_x(0,0) = 1;
    f_x(0,3) = dt;
    f_x(1,1) = 1;
    f_x(1,4) = dt;
    f_x(2,2) = 1;
    f_x(2,5) = dt;
    f_x(3,3) = 1;
    f_x(3,6) = -x11*(-x2 + x3*x6);
    f_x(3,7) = x12*x14;
    f_x(3,8) = x11*(-x2*x3 + x6);
    f_x(4,4) = 1;
    f_x(4,6) = x11*(x12*x3 + x15);
    f_x(4,7) = x14*x6;
    f_x(4,8) = -x11*(x12 + x15*x3);
    f_x(5,5) = 1;
    f_x(5,7) = -x11*x3*x5;
    f_x(5,8) = -x0*x14;
    f_x(6,6) = 1;
    f_x(6,9) = dt;
    f_x(7,7) = 1;
    f_x(7,10) = dt;
    f_x(8,8) = 1;
    f_x(8,11) = dt;
    f_x(9,6) = dt*(x16*(-this->m_I_yy*x20 + this->m_I_yy*x24 + this->m_I_zz*x20 - this->m_I_zz*x24) + x35*x41 - x41*x49 + x41*x55 - x41*x59);
    f_x(9,7) = dt*(-x1*x71*x72 - x16*x25*(x61 + x63 - x64 - x65) - x25*x[10] + x4*x72*x76 + x78*x79 + x78 + x79*x81 + x81);
    f_x(9,9) = x3*x87 + 1;
    f_x(9,10) = dt*(-x16*(this->m_I_yy*x92 + this->m_I_yy*x93 - this->m_I_zz*x92 - this->m_I_zz*x93) + x18 + x41*x90 + x41*x91);
    f_x(9,11) = dt*(x101*x41 - x13*x16*(-x61 - x63 + x64 + x65) + x13*x[10] + x41*x98);
    f_x(10,6) = -x85*(-x1*(-x19*x27*x68*(-this->m_I_xx + this->m_I_zz) + x34) - x4*(-x43*x68*x83 + x48) + x77 + x80);
    f_x(10,7) = x85*(x102*x71 + x103*x76);
    f_x(10,9) = -x85*(x1*x84 + x4*x82);
    f_x(10,10) = x1*x4*x85*(x44 + x88) + 1;
    f_x(10,11) = x85*(x1*x100 - x4*x97);
    f_x(11,6) = x86*(x35 - x49 + x55 - x59);
    f_x(11,7) = x86*(x102*x76 - x103*x71 + x40*x77 + x40*x80);
    f_x(11,9) = x87;
    f_x(11,10) = x86*(x90 + x91);
    f_x(11,11) = x86*(x101 + x98) + 1;

    // Return the Jacobian
    return f_x;
}

arma::mat QuadcopterDynamics::f_u(arma::vec x, arma::vec u, double t, double dt) {
    // Common subexpressions
    double x0 = sqrt(this->m_g*this->m_m/this->m_k);
    double x1 = dt*(2*u[0] + x0);
    double x2 = sin(x[6]);
    double x3 = sin(x[8]);
    double x4 = cos(x[6]);
    double x5 = cos(x[8]);
    double x6 = sin(x[7]);
    double x7 = x5*x6;
    double x8 = this->m_k/this->m_m;
    double x9 = x8*(x2*x3 + x4*x7);
    double x10 = 2*u[1] + x0;
    double x11 = dt*x9;
    double x12 = dt*(2*u[2] + x0);
    double x13 = 2*u[3] + x0;
    double x14 = x8*(x2*x7 - x3*x4);
    double x15 = dt*x14;
    double x16 = cos(x[7]);
    double x17 = x16*x5*x8;
    double x18 = dt*x17;
    double x19 = this->m_k*this->m_l;
    double x20 = x19/this->m_I_xx;
    double x21 = 1.0/(pow(x2, 2) + pow(x4, 2));
    double x22 = 1.0/x16;
    double x23 = this->m_b/this->m_I_zz;
    double x24 = x23*x4;
    double x25 = x21*x22*x24;
    double x26 = x25*x6;
    double x27 = x19/this->m_I_yy;
    double x28 = x2*x27;
    double x29 = dt*x21;
    double x30 = x10*x29;
    double x31 = x22*x30*(x24 - x28);
    double x32 = x13*x29;
    double x33 = x22*x32*(x24 + x28);
    double x34 = x2*x23;
    double x35 = x21*x34;
    double x36 = x27*x4;

    // Allocate matrix for Jacobian with respect to the control
    arma::mat f_u = arma::zeros<arma::mat>(x.n_elem, u.n_elem);

    // Compute components of Jacobian
    f_u(3,0) = x1*x9;
    f_u(3,1) = x10*x11;
    f_u(3,2) = x12*x9;
    f_u(3,3) = x11*x13;
    f_u(4,0) = x1*x14;
    f_u(4,1) = x10*x15;
    f_u(4,2) = x12*x14;
    f_u(4,3) = x13*x15;
    f_u(5,0) = x1*x17;
    f_u(5,1) = x10*x18;
    f_u(5,2) = x12*x17;
    f_u(5,3) = x13*x18;
    f_u(9,0) = x1*(x20 + x26);
    f_u(9,1) = -x31*x6;
    f_u(9,2) = x12*(-x20 + x26);
    f_u(9,3) = -x33*x6;
    f_u(10,0) = -x1*x35;
    f_u(10,1) = x30*(x34 + x36);
    f_u(10,2) = -x12*x35;
    f_u(10,3) = -x32*(-x34 + x36);
    f_u(11,0) = x1*x25;
    f_u(11,1) = -x31;
    f_u(11,2) = x12*x25;
    f_u(11,3) = -x33;

    // Return the Jacobian
    return f_u;
}

arma::cube QuadcopterDynamics::f_xx(arma::vec x, arma::vec u, double t, double dt) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = x0*x1;
    double x3 = sin(x[7]);
    double x4 = cos(x[6]);
    double x5 = cos(x[8]);
    double x6 = x4*x5;
    double x7 = x3*x6;
    double x8 = sqrt(this->m_g*this->m_m/this->m_k);
    double x9 = pow(2*u[1] + x8, 2);
    double x10 = pow(2*u[3] + x8, 2);
    double x11 = pow(2*u[0] + x8, 2) + pow(2*u[2] + x8, 2);
    double x12 = (1.0/4.0)*dt*this->m_k*(x10 + x11 + x9)/this->m_m;
    double x13 = -x12*(x2 + x7);
    double x14 = x0*x5;
    double x15 = cos(x[7]);
    double x16 = x12*x15;
    double x17 = -x14*x16;
    double x18 = x12*(x2*x3 + x6);
    double x19 = x1*x4;
    double x20 = -x16*x19;
    double x21 = x14*x3;
    double x22 = -x12*(-x19 + x21);
    double x23 = x16*x6;
    double x24 = x12*(x14 - x19*x3);
    double x25 = -x16*x2;
    double x26 = -x16*x5;
    double x27 = x1*x12*x3;
    double x28 = x4*x[10];
    double x29 = x15*x[11];
    double x30 = x0*x29;
    double x31 = x28 + x30;
    double x32 = x0*x[10];
    double x33 = x29*x4;
    double x34 = x32 - x33;
    double x35 = 1.0/this->m_I_xx;
    double x36 = -this->m_I_zz;
    double x37 = this->m_I_yy + x36;
    double x38 = x35*x37;
    double x39 = x3*x[11];
    double x40 = -x39;
    double x41 = x40 + x[9];
    double x42 = 1.0/this->m_I_zz;
    double x43 = -this->m_I_xx;
    double x44 = x42*(this->m_I_yy + x43);
    double x45 = x41*x44;
    double x46 = x15*x[9];
    double x47 = x0*x46 + x28*x3;
    double x48 = x47*x[11];
    double x49 = x28*x[9];
    double x50 = x48 + x49;
    double x51 = -x31*x45 + x50;
    double x52 = x4*x51;
    double x53 = pow(x0, 2);
    double x54 = pow(x4, 2);
    double x55 = 1.0/(x53 + x54);
    double x56 = 1.0/x15;
    double x57 = x3*x56;
    double x58 = x55*x57;
    double x59 = 1.0/this->m_I_yy;
    double x60 = x59*(this->m_I_zz + x43);
    double x61 = x31*x60;
    double x62 = -x41*x61 + x50;
    double x63 = x4*x55;
    double x64 = 2*x57;
    double x65 = x3*x32;
    double x66 = x4*x46;
    double x67 = x65 - x66;
    double x68 = x67*x[11];
    double x69 = x32*x[9];
    double x70 = x68 + x69;
    double x71 = -x34*x41*x60 + x70;
    double x72 = x0*x71;
    double x73 = -x32 + x33;
    double x74 = x42*(this->m_I_xx - this->m_I_yy);
    double x75 = x41*x74;
    double x76 = x70 - x73*x75;
    double x77 = x0*x55;
    double x78 = x39 - x[9];
    double x79 = 4*x78;
    double x80 = x34*x79;
    double x81 = -x10;
    double x82 = x59*(-this->m_I_xx*x80 + this->m_I_zz*x80 + this->m_k*this->m_l*(x81 + x9)) + 4*x68 + 4*x69;
    double x83 = (1.0/4.0)*x82;
    double x84 = x57*x77;
    double x85 = x31*x79;
    double x86 = x42*(-this->m_I_xx*x85 + this->m_I_yy*x85 + this->m_b*(x11 + x81 - x9)) + 4*x48 + 4*x49;
    double x87 = (1.0/4.0)*x86;
    double x88 = x57*x63;
    double x89 = x62*x77;
    double x90 = x63*x76;
    double x91 = pow(x15, 2);
    double x92 = 1.0/x91;
    double x93 = pow(x3, 2);
    double x94 = x92*x93;
    double x95 = x0*x31;
    double x96 = this->m_I_zz*x95;
    double x97 = x4*x73;
    double x98 = this->m_I_yy*x95;
    double x99 = 2*x35;
    double x100 = x15*x31;
    double x101 = this->m_I_zz*x100;
    double x102 = x3*x41;
    double x103 = this->m_I_zz*x0;
    double x104 = this->m_I_xx*x100;
    double x105 = -x104;
    double x106 = this->m_I_xx*x102;
    double x107 = x0*x106;
    double x108 = x105 - x107;
    double x109 = x3*x[9];
    double x110 = x0*x109;
    double x111 = x15*x28;
    double x112 = x110 - x111;
    double x113 = x112 - x59*(x101 + x102*x103 + x108);
    double x114 = x39*x56;
    double x115 = x114*x77;
    double x116 = this->m_I_yy*x102;
    double x117 = x116*x4;
    double x118 = x106*x4;
    double x119 = x15*x73;
    double x120 = this->m_I_yy*x119;
    double x121 = this->m_I_xx*x119;
    double x122 = x109*x4 + x15*x32;
    double x123 = x122 - x42*(x117 - x118 + x120 - x121);
    double x124 = x114*x63;
    double x125 = this->m_I_yy*x100;
    double x126 = x0*x116;
    double x127 = x112 - x42*(x108 + x125 + x126);
    double x128 = x4*x78;
    double x129 = x128*x3;
    double x130 = this->m_I_zz*x129;
    double x131 = this->m_I_xx*x129;
    double x132 = x15*x34;
    double x133 = this->m_I_zz*x132;
    double x134 = this->m_I_xx*x132;
    double x135 = x122 + x59*(x130 - x131 + x133 - x134);
    double x136 = x63*x83;
    double x137 = x77*x87;
    double x138 = -dt*(x113*x115 - x115*x127 + x123*x124 - x124*x135 - x136*x94 - x136 + x137*x94 + x137 + x39*x99*(this->m_I_yy*x97 - this->m_I_zz*x97 + x96 - x98) - x89*x94 - x89 + x90*x94 + x90);
    double x139 = x31 - x61;
    double x140 = x31*x74 + x31;
    double x141 = x59*(this->m_I_xx + x36);
    double x142 = x141*x73;
    double x143 = x142 + x73;
    double x144 = x34*x74;
    double x145 = x0*x139 - x0*x140 - x143*x4 + x4*(-x144 + x73);
    double x146 = dt*x55;
    double x147 = x146*x3;
    double x148 = x147*x56;
    double x149 = x145*x148;
    double x150 = x39 + x[9];
    double x151 = x150 - x45;
    double x152 = x0*x57;
    double x153 = x152*x63;
    double x154 = x150 + x44*x78;
    double x155 = x150 + x60*x78;
    double x156 = 2*x0;
    double x157 = x0*x73;
    double x158 = x31*x4;
    double x159 = this->m_I_yy*x158 - this->m_I_zz*x158;
    double x160 = -dt*(x151*x153 + x153*x154 - x155*x156*x88 + x99*(this->m_I_yy*x157 - this->m_I_zz*x157 + x159));
    double x161 = x34*x4;
    double x162 = x3*x31;
    double x163 = this->m_I_zz*x162;
    double x164 = x15*x78;
    double x165 = x103*x164;
    double x166 = this->m_I_xx*x162;
    double x167 = x0*x164;
    double x168 = this->m_I_xx*x167;
    double x169 = x166 + x168;
    double x170 = x0*(x47 - x59*(-x163 - x165 + x169));
    double x171 = x3*x73;
    double x172 = x15*x41;
    double x173 = x172*x4;
    double x174 = this->m_I_xx*x171;
    double x175 = this->m_I_xx*x173;
    double x176 = x174 - x175;
    double x177 = x42*(-this->m_I_yy*x171 + this->m_I_yy*x173 + x176);
    double x178 = x4*(x177 + x67);
    double x179 = -this->m_I_yy*x162;
    double x180 = -this->m_I_yy*x167 + x169 + x179;
    double x181 = -x180*x42 + x47;
    double x182 = x0*x181;
    double x183 = x128*x15;
    double x184 = this->m_I_zz*x183;
    double x185 = this->m_I_xx*x183;
    double x186 = x3*x34;
    double x187 = this->m_I_xx*x186;
    double x188 = this->m_I_zz*x186;
    double x189 = -x59*(x184 - x185 + x187 - x188) + x67;
    double x190 = x189*x4;
    double x191 = dt*(-x15*x99*(this->m_I_yy*x161 - this->m_I_zz*x161 - x96 + x98) + x170*x58 - x178*x58 - x182*x58 + x190*x58);
    double x192 = x112 + x42*(x104 + x107 - x125 - x126);
    double x193 = 2*x[11];
    double x194 = x192*x193*x63;
    double x195 = -x130 + x131 - x133 + x134;
    double x196 = x122 - x195*x59;
    double x197 = x193*x196*x77;
    double x198 = 2*x4;
    double x199 = x0*x[11];
    double x200 = this->m_I_yy*x199;
    double x201 = x200*x93;
    double x202 = this->m_I_zz*x198;
    double x203 = this->m_I_yy*x0;
    double x204 = this->m_I_zz*x171;
    double x205 = this->m_I_zz*x173;
    double x206 = x15*x39;
    double x207 = x202*x206;
    double x208 = this->m_I_xx*x206;
    double x209 = x198*x208;
    double x210 = x59*(x176 - x204 + x205 - x207 + x209) + x67;
    double x211 = -this->m_I_yy*x156*x206 + x156*x208;
    double x212 = -x42*(x180 + x211) + x47;
    double x213 = (1.0/2.0)*x82;
    double x214 = pow(x3, 3)/pow(x15, 3);
    double x215 = (1.0/2.0)*x86;
    double x216 = x141 + 1;
    double x217 = x199*x4;
    double x218 = x216*x217;
    double x219 = x56*x93;
    double x220 = x74 + 1;
    double x221 = x217*x220;
    double x222 = x140*x4;
    double x223 = x0*x143;
    double x224 = x146*(x218*x219 - x219*x221 + x222*x94 + x222 - x223*x94 - x223);
    double x225 = -x141;
    double x226 = x225 + 1;
    double x227 = x226*x53;
    double x228 = x39*x55;
    double x229 = -x74;
    double x230 = x229 + 1;
    double x231 = x230*x54;
    double x232 = x35*(this->m_I_yy*x53 - this->m_I_yy*x54 - this->m_I_zz*x53 + this->m_I_zz*x54);
    double x233 = x141*x78;
    double x234 = x53*(x150 - x233);
    double x235 = x234*x55;
    double x236 = x74*x78;
    double x237 = x54*(x150 - x236);
    double x238 = x237*x55;
    double x239 = dt*(x227*x228 + x228*x231 + x232*x39 + x235*x94 + x235 + x238*x94 + x238 + x40);
    double x240 = x0*x34;
    double x241 = x181*x63;
    double x242 = x189*x77;
    double x243 = this->m_I_xx*x199;
    double x244 = x3*x78;
    double x245 = this->m_I_xx*x0;
    double x246 = x105 + x244*x245;
    double x247 = -x110 + x111;
    double x248 = x247 + x42*(x125 + x200*x91 - x201 - x203*x244 - x243*x91 + x243*x93 + x246);
    double x249 = x248*x4;
    double x250 = x4*x[11];
    double x251 = this->m_I_xx*x250;
    double x252 = this->m_I_zz*x250;
    double x253 = x122 - x59*(x195 - x251*x91 + x251*x93 + x252*x91 - x252*x93);
    double x254 = x0*x253;
    double x255 = dt*(x241*x94 + x241 + x242*x94 + x242 + x249*x58 + x254*x58 - x3*x35*(this->m_I_yy*x198*x30 - this->m_I_yy*x240 + this->m_I_zz*x240 + x159 - x202*x30) - x3*x[10]);
    double x256 = x216*x53 + x220*x54;
    double x257 = x148*x256;
    double x258 = x225 + x74;
    double x259 = dt*x0;
    double x260 = x259*x63;
    double x261 = x258*x260;
    double x262 = x261*x3;
    double x263 = x55*x93;
    double x264 = x263*x56;
    double x265 = dt*(-x15*x232 + x15 + x227*x264 + x231*x264);
    double x266 = (1.0/4.0)*x4;
    double x267 = x266*x82;
    double x268 = (1.0/4.0)*x0;
    double x269 = x268*x86;
    double x270 = -x146*(x113*x250 - x123*x199 - x127*x250 + x135*x199);
    double x271 = -x146*(x0*(-x142 + x34) + x0*(x34*x44 + x73) - x139*x4 + x222);
    double x272 = -x146*(-x151*x53 + x154*x54 + x155*x53 - x155*x54);
    double x273 = x166 + x172*x203 - x172*x245 + x179;
    double x274 = -x146*(x0*(x59*(-x184 + x185 - x187 + x188) + x67) + x0*(-x177 - x65 + x66) - x4*(x47 + x59*(x163 + x165 - x166 - x168)) + x4*(-x273*x42 + x47));
    double x275 = x220*x53;
    double x276 = x147*(x216*x54*x[11] + x275*x[11]);
    double x277 = x260*(x226*x29 - x230*x29);
    double x278 = x146*(-x0*x248 + x253*x4);
    double x279 = x260*(x141 + x229);
    double x280 = -x146*x15*(x275 + x54*(-x60 + 1));
    double x281 = x266*x86 + x268*x82;
    double x282 = x146*x56;
    double x283 = x4*x57;
    double x284 = x282*(x152*(-x233*x31 + x50) + x192*x199 + x196*x250 + x199*(x247 + x59*(x101 - x103*x244 + x246)) - x250*(x122 + x42*(-x117 + x118 - x120 + x121)) + x267*x57 - x269*x57 - x283*(-x144*x78 + x70));
    double x285 = x145*x282;
    double x286 = x260*x56*(-2*x233 + x236 - x75);
    double x287 = x282*(x170 - x178 - x182 + x190);
    double x288 = (1.0/2.0)*x94;
    double x289 = x148*(x218 - x221 + x222*x56 - x223*x56);
    double x290 = x3*x92;
    double x291 = x146*(x227*x[11] + x231*x[11] + x234*x290 + x237*x290);
    double x292 = x282*(x152*x189 + x181*x283 + x249 + x254);
    double x293 = x256*x282;
    double x294 = x148*(x227 + x231);

    // Allocate cube for second derivative of f with respect to the state twice
    arma::cube f_xx = arma::zeros<arma::cube>(x.n_elem, x.n_elem, x.n_elem);

    // Compute components of tensor
    f_xx(3,6,6) = x13;
    f_xx(3,6,7) = x17;
    f_xx(3,6,8) = x18;
    f_xx(3,7,6) = x17;
    f_xx(3,7,7) = -x12*x7;
    f_xx(3,7,8) = x20;
    f_xx(3,8,6) = x18;
    f_xx(3,8,7) = x20;
    f_xx(3,8,8) = x13;
    f_xx(4,6,6) = x22;
    f_xx(4,6,7) = x23;
    f_xx(4,6,8) = x24;
    f_xx(4,7,6) = x23;
    f_xx(4,7,7) = -x12*x21;
    f_xx(4,7,8) = x25;
    f_xx(4,8,6) = x24;
    f_xx(4,8,7) = x25;
    f_xx(4,8,8) = x22;
    f_xx(5,7,7) = x26;
    f_xx(5,7,8) = x27;
    f_xx(5,8,7) = x27;
    f_xx(5,8,8) = x26;
    f_xx(9,6,6) = -dt*(-4*x31*x34*x38 + x52*x58 + x58*x72 - x62*x63*x64 - x64*x76*x77 + x83*x84 + x87*x88);
    f_xx(9,6,7) = x138;
    f_xx(9,6,9) = x149;
    f_xx(9,6,10) = x160;
    f_xx(9,6,11) = x191;
    f_xx(9,7,6) = x138;
    f_xx(9,7,7) = dt*(-x115*x210 - x124*x212 - x194*x94 - x194 + x197*x94 + x197 + x213*x214*x77 + x213*x84 + x214*x215*x63 + x215*x88 - x29*x[10] + x35*x[11]*(-x0*x133 + x101*x4 - x125*x4 + x132*x203 + x198*x201 - x199*x202*x93));
    f_xx(9,7,9) = x224;
    f_xx(9,7,10) = x239;
    f_xx(9,7,11) = x255;
    f_xx(9,9,6) = x149;
    f_xx(9,9,7) = x224;
    f_xx(9,9,10) = x257;
    f_xx(9,9,11) = x262;
    f_xx(9,10,6) = x160;
    f_xx(9,10,7) = x239;
    f_xx(9,10,9) = x257;
    f_xx(9,10,10) = -x259*x37*x4*x99;
    f_xx(9,10,11) = x265;
    f_xx(9,11,6) = x191;
    f_xx(9,11,7) = x255;
    f_xx(9,11,9) = x262;
    f_xx(9,11,10) = x265;
    f_xx(9,11,11) = x198*x259*(x141*x263 - x263*x74 + x38*x91);
    f_xx(10,6,6) = -x146*(-x0*x51 + x156*x62 - x198*x76 + x267 - x269 + x4*x71);
    f_xx(10,6,7) = x270;
    f_xx(10,6,9) = x271;
    f_xx(10,6,10) = x272;
    f_xx(10,6,11) = x274;
    f_xx(10,7,6) = x270;
    f_xx(10,7,7) = -x146*(-x199*(-x42*(x211 + x273) + x47) + x250*(-x59*(-x174 + x175 + x204 - x205 + x207 - x209) + x67));
    f_xx(10,7,9) = x276;
    f_xx(10,7,10) = x277;
    f_xx(10,7,11) = x278;
    f_xx(10,9,6) = x271;
    f_xx(10,9,7) = x276;
    f_xx(10,9,10) = x279;
    f_xx(10,9,11) = x280;
    f_xx(10,10,6) = x272;
    f_xx(10,10,7) = x277;
    f_xx(10,10,9) = x279;
    f_xx(10,10,11) = x262;
    f_xx(10,11,6) = x274;
    f_xx(10,11,7) = x278;
    f_xx(10,11,9) = x280;
    f_xx(10,11,10) = x262;
    f_xx(10,11,11) = 2*x147*x15*(x141*x54 + x53*x74);
    f_xx(11,6,6) = -x282*(-x156*x76 - x198*x62 + x281 + x52 + x72);
    f_xx(11,6,7) = x284;
    f_xx(11,6,9) = x285;
    f_xx(11,6,10) = x286;
    f_xx(11,6,11) = x287;
    f_xx(11,7,6) = x284;
    f_xx(11,7,7) = x282*(x0*x288*x82 + x114*x156*x196 - x114*x192*x198 - x199*x210 - x212*x250 + x281 + x288*x4*x86);
    f_xx(11,7,9) = x289;
    f_xx(11,7,10) = x291;
    f_xx(11,7,11) = x292;
    f_xx(11,9,6) = x285;
    f_xx(11,9,7) = x289;
    f_xx(11,9,10) = x293;
    f_xx(11,9,11) = x261;
    f_xx(11,10,6) = x286;
    f_xx(11,10,7) = x291;
    f_xx(11,10,9) = x293;
    f_xx(11,10,11) = x294;
    f_xx(11,11,6) = x287;
    f_xx(11,11,7) = x292;
    f_xx(11,11,9) = x261;
    f_xx(11,11,10) = x294;
    f_xx(11,11,11) = -dt*x156*x258*x3*x63;

    // Return the tensor
    return f_xx;
}

arma::cube QuadcopterDynamics::f_uu(arma::vec x, arma::vec u, double t, double dt) {
    // Common subexpressions
    double x0 = sin(x[6]);
    double x1 = sin(x[8]);
    double x2 = cos(x[6]);
    double x3 = cos(x[8]);
    double x4 = sin(x[7]);
    double x5 = x3*x4;
    double x6 = 2*dt;
    double x7 = this->m_k*x6/this->m_m;
    double x8 = x7*(x0*x1 + x2*x5);
    double x9 = x7*(x0*x5 - x1*x2);
    double x10 = cos(x[7]);
    double x11 = x10*x3*x7;
    double x12 = this->m_k*this->m_l;
    double x13 = x12/this->m_I_xx;
    double x14 = 1.0/(pow(x0, 2) + pow(x2, 2));
    double x15 = 1.0/x10;
    double x16 = this->m_b/this->m_I_zz;
    double x17 = x16*x2;
    double x18 = x14*x15*x17*x4;
    double x19 = x12/this->m_I_yy;
    double x20 = x0*x19;
    double x21 = x14*x6;
    double x22 = x15*x21;
    double x23 = x22*(x17 - x20);
    double x24 = x22*(x17 + x20);
    double x25 = x0*x16;
    double x26 = -x21*x25;
    double x27 = x19*x2;
    double x28 = x17*x22;

    // Allocate cube for second derivative of f with respect to the control twice
    arma::cube f_uu = arma::zeros<arma::cube>(x.n_elem, u.n_elem, u.n_elem);

    // Compute components of tensor
    f_uu(3,0,0) = x8;
    f_uu(3,1,1) = x8;
    f_uu(3,2,2) = x8;
    f_uu(3,3,3) = x8;
    f_uu(4,0,0) = x9;
    f_uu(4,1,1) = x9;
    f_uu(4,2,2) = x9;
    f_uu(4,3,3) = x9;
    f_uu(5,0,0) = x11;
    f_uu(5,1,1) = x11;
    f_uu(5,2,2) = x11;
    f_uu(5,3,3) = x11;
    f_uu(9,0,0) = x6*(x13 + x18);
    f_uu(9,1,1) = -x23*x4;
    f_uu(9,2,2) = x6*(-x13 + x18);
    f_uu(9,3,3) = -x24*x4;
    f_uu(10,0,0) = x26;
    f_uu(10,1,1) = x21*(x25 + x27);
    f_uu(10,2,2) = x26;
    f_uu(10,3,3) = x21*(x25 - x27);
    f_uu(11,0,0) = x28;
    f_uu(11,1,1) = -x23;
    f_uu(11,2,2) = x28;
    f_uu(11,3,3) = -x24;

    // Return the tensor
    return f_uu;
}

arma::cube QuadcopterDynamics::f_ux(arma::vec x, arma::vec u, double t, double dt) {
    // Common subexpressions
    double x0 = sin(x[8]);
    double x1 = cos(x[6]);
    double x2 = x0*x1;
    double x3 = sin(x[7]);
    double x4 = sin(x[6]);
    double x5 = cos(x[8]);
    double x6 = x4*x5;
    double x7 = -x2 + x3*x6;
    double x8 = sqrt(this->m_g*this->m_m/this->m_k);
    double x9 = 2*u[0] + x8;
    double x10 = dt*this->m_k/this->m_m;
    double x11 = x10*x9;
    double x12 = cos(x[7]);
    double x13 = x1*x5;
    double x14 = x12*x13;
    double x15 = -x2*x3 + x6;
    double x16 = 2*u[1] + x8;
    double x17 = x10*x16;
    double x18 = 2*u[2] + x8;
    double x19 = x10*x18;
    double x20 = 2*u[3] + x8;
    double x21 = x10*x20;
    double x22 = x0*x4;
    double x23 = x13*x3 + x22;
    double x24 = x12*x6;
    double x25 = x13 + x22*x3;
    double x26 = x3*x5;
    double x27 = x0*x12;
    double x28 = 1.0/x12;
    double x29 = this->m_b/this->m_I_zz;
    double x30 = x29*x4;
    double x31 = dt/(pow(x1, 2) + pow(x4, 2));
    double x32 = x28*x30*x31;
    double x33 = x32*x9;
    double x34 = pow(x12, -2);
    double x35 = pow(x3, 2)*x34;
    double x36 = x35 + 1;
    double x37 = x1*x29;
    double x38 = x31*x37;
    double x39 = x38*x9;
    double x40 = this->m_k*this->m_l/this->m_I_yy;
    double x41 = x1*x40;
    double x42 = x16*x31;
    double x43 = x28*x42*(x30 + x41);
    double x44 = x35*x37;
    double x45 = x4*x40;
    double x46 = x35*x45;
    double x47 = x37 - x45;
    double x48 = x18*x32;
    double x49 = x18*x38;
    double x50 = x20*x31;
    double x51 = x28*x50*(-x30 + x41);
    double x52 = x37 + x45;
    double x53 = x50*x52;
    double x54 = x3*x34;

    // Allocate cube for second derivative of f with respect to the control and then state
    arma::cube f_ux = arma::zeros<arma::cube>(x.n_elem, u.n_elem, x.n_elem);

    // Compute components of tensor
    f_ux(3,0,6) = -x11*x7;
    f_ux(3,0,7) = x11*x14;
    f_ux(3,0,8) = x11*x15;
    f_ux(3,1,6) = -x17*x7;
    f_ux(3,1,7) = x14*x17;
    f_ux(3,1,8) = x15*x17;
    f_ux(3,2,6) = -x19*x7;
    f_ux(3,2,7) = x14*x19;
    f_ux(3,2,8) = x15*x19;
    f_ux(3,3,6) = -x21*x7;
    f_ux(3,3,7) = x14*x21;
    f_ux(3,3,8) = x15*x21;
    f_ux(4,0,6) = x11*x23;
    f_ux(4,0,7) = x11*x24;
    f_ux(4,0,8) = -x11*x25;
    f_ux(4,1,6) = x17*x23;
    f_ux(4,1,7) = x17*x24;
    f_ux(4,1,8) = -x17*x25;
    f_ux(4,2,6) = x19*x23;
    f_ux(4,2,7) = x19*x24;
    f_ux(4,2,8) = -x19*x25;
    f_ux(4,3,6) = x21*x23;
    f_ux(4,3,7) = x21*x24;
    f_ux(4,3,8) = -x21*x25;
    f_ux(5,0,7) = -x11*x26;
    f_ux(5,0,8) = -x11*x27;
    f_ux(5,1,7) = -x17*x26;
    f_ux(5,1,8) = -x17*x27;
    f_ux(5,2,7) = -x19*x26;
    f_ux(5,2,8) = -x19*x27;
    f_ux(5,3,7) = -x21*x26;
    f_ux(5,3,8) = -x21*x27;
    f_ux(9,0,6) = -x3*x33;
    f_ux(9,0,7) = x36*x39;
    f_ux(9,1,6) = x3*x43;
    f_ux(9,1,7) = -x42*(x44 - x46 + x47);
    f_ux(9,2,6) = -x3*x48;
    f_ux(9,2,7) = x36*x49;
    f_ux(9,3,6) = -x3*x51;
    f_ux(9,3,7) = -x50*(x44 + x46 + x52);
    f_ux(10,0,6) = -x39;
    f_ux(10,1,6) = -x42*(-x37 + x45);
    f_ux(10,2,6) = -x49;
    f_ux(10,3,6) = x53;
    f_ux(11,0,6) = -x33;
    f_ux(11,0,7) = x39*x54;
    f_ux(11,1,6) = x43;
    f_ux(11,1,7) = -x42*x47*x54;
    f_ux(11,2,6) = -x48;
    f_ux(11,2,7) = x49*x54;
    f_ux(11,3,6) = -x51;
    f_ux(11,3,7) = -x53*x54;

    // Return the tensor
    return f_ux;
}