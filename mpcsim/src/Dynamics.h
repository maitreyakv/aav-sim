/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef MPCSIM_DYNAMICS_H
#define MPCSIM_DYNAMICS_H


#include <armadillo>

// TODO: Update documentation

class Dynamics {

public:

    virtual arma::vec F(arma::vec x, arma::vec u, double t) = 0;

    virtual arma::vec f(arma::vec x, arma::vec u, double t, double dt) = 0;

    virtual arma::mat f_x(arma::vec x, arma::vec u, double t, double dt) = 0;

    virtual arma::mat f_u(arma::vec x, arma::vec u, double t, double dt) = 0;

    virtual arma::cube f_xx(arma::vec x, arma::vec u, double t, double dt) = 0;

    virtual arma::cube f_uu(arma::vec x, arma::vec u, double t, double dt) = 0;

    virtual arma::cube f_ux(arma::vec x, arma::vec u, double t, double dt) = 0;

    /**
     * Gets the dimensionality (size) of the state vector for the Dynamics
     *
     * @return The size of the state vector
     */
    virtual int getStateDimension() = 0;

    /**
     * Gets the dimensionality (size) of the control vector for the Dynamics
     *
     * @return The size of the control vector
     */
    virtual int getControlDimension() = 0;
};


#endif //MPCSIM_DYNAMICS_H
