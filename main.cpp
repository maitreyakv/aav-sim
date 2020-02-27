#include <iostream>
#include "src/CartPoleDynamics.h"
#include "src/DDPController.h"
#include "src/QuadraticCost.h"
#include "src/Simulation.h"

const double pi = atan(1.0) * 4.0;

int main() {

    CartPoleDynamics* cart_pole_dynamics_ptr = new CartPoleDynamics(1.0, 0.01, 0.25);

    arma::mat Q_f = 10.0 * arma::eye<arma::mat>(4, 4);
    arma::mat R = 0.1 * arma::eye<arma::mat>(1, 1);
    QuadraticCost* quadratic_cost_ptr = new QuadraticCost(Q_f, R);

    DDPController* ddp_controller_ptr = new DDPController(cart_pole_dynamics_ptr,quadratic_cost_ptr, 200, 100, 0.1);

    arma::vec x_0;
    x_0 << 0.0 << arma::endr << 0.0 << arma::endr << 0.1 << arma::endr << 0.0 << arma::endr;

    arma::vec x_star;
    x_star << 0.0 << arma::endr << 0.0 << arma::endr << pi << arma::endr << 0.0 << arma::endr;

    System* system_ptr = new System(cart_pole_dynamics_ptr, ddp_controller_ptr);

    Simulation* simulation_ptr = new Simulation(system_ptr, 0.05, x_0, x_star, 0.5);

    simulation_ptr->simulate(3.0);

    return 0;
}