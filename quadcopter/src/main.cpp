#include <iostream>
#include "QuadcopterDynamics.h"
#include "../../mpcsim/src/DDPController.h"
#include "../../mpcsim/src/QuadraticCost.h"
#include "../../mpcsim/src/Simulation.h"

const double pi = atan(1.0) * 4.0;

int main() {

    double m = 0.468;
    double g = 9.81;
    double l = 0.225;
    double k = 0.00000298;
    double b = 0.000000114;
    double I_xx = 0.004856;
    double I_yy = 0.004856;
    double I_zz = 0.008801;

    QuadcopterDynamics* quadcopter_dynamics_ptr = new QuadcopterDynamics(m, g, l, k, b, I_xx, I_yy, I_zz);

    arma::vec x = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    arma::vec x_star = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

/**
    double dt = 0.01;

    arma::vec u = {1.0, 1.0, 1.0, 1.0};

    for (int k = 0; k < 1000; k++) {
        x = x + quadcopter_dynamics_ptr->F(x, u, k*dt) * dt;
        std::cout << k*dt << "," << u[0] << "," << u[1] << "," << u[2] << "," << u[3];
        for (int i = 0; i < 12; i++) {
            std::cout << "," << x[i];
        }
        std::cout << std::endl;
    }
*/

    arma::mat Q_f = 10.0 * arma::eye<arma::mat>(12, 12);
    arma::mat R = 0.01 * arma::eye<arma::mat>(4, 4);
    QuadraticCost* quadratic_cost_ptr = new QuadraticCost(Q_f, R);

    arma::vec u_max = pow(10, 1) * arma::ones<arma::vec>(4);

    DDPController* ddp_controller_ptr = new DDPController(quadcopter_dynamics_ptr, quadratic_cost_ptr, u_max, 100, 50, 0.5);

    System* system_ptr = new System(quadcopter_dynamics_ptr, ddp_controller_ptr);

    Simulation* simulation_ptr = new Simulation(system_ptr, 0.02, x, x_star, 2.0);

    simulation_ptr->simulate(5.0);

    return 0;
}