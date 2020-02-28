#include <iostream>
#include "src/QuadcopterDynamics.h"
#include "src/DDPController.h"
#include "src/QuadraticCost.h"
#include "src/Simulation.h"

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

    double dt = 0.01;

    double u_hover = 1.01*(m*g/4)/k;
    arma::vec u = {1.001*u_hover, u_hover, u_hover, u_hover};

    for (int k = 0; k < 1000; k++) {
        x = x + quadcopter_dynamics_ptr->F(x, u, k*dt) * dt;
        std::cout << k*dt << "," << u[0] << "," << u[1] << "," << u[2] << "," << u[3];
        for (int i = 0; i < 12; i++) {
            std::cout << "," << x[i];
        }
        std::cout << std::endl;
    }

    return 0;
}