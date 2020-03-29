/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <fstream>
#include <boost/numeric/odeint.hpp>

#include "Simulation.h"
#include "../navigation/RRTPathPlanner.h"

// Required for using Armadillo Vectors with odeint
namespace boost {
    namespace numeric {
        namespace odeint {
            template<>
            struct is_resizeable< arma::vec >
            {
                typedef boost::true_type type;
                static const bool value = type::value;
            };
        }
    }
}

int Simulation::simulate(double time) {

    // TEMP
    std::vector<Obstacle*> obs;
    RRTPathPlanner* path_planner_ptr = new RRTPathPlanner(obs);
    arma::vec pos_start = {this->m_x(0), this->m_x(1), this->m_x(2)};
    arma::vec pos_target = {this->m_x_star(0), this->m_x_star(1), this->m_x_star(2)};
    std::vector<arma::vec> path = path_planner_ptr->computePath(pos_start, pos_target);
    for (int i = 0; i < path.size(); i++) {
        std::cout << path[i](0) << "," << path[i](1) << "," << path[i](2) << std::endl;
    }
    delete path_planner_ptr;


    // Make armadillo vectors resizable in odeint
    BOOST_STATIC_ASSERT(boost::numeric::odeint::is_resizeable<arma::vec>::value == true);

    // Initialize the time horizon with the specified horizon, this may be edited during the simulation
    double horizon = this->m_horizon;

    // Open a text file for the simulation output
    std::ofstream output_file;
    output_file.open("output.txt", std::ios::trunc);

    // Perform the simulation by updating the control and integrating the system
    while (this->m_t < time) {
        // Obtain the current control
        arma::vec u = this->m_system_ptr->getControl();

        // Write the current time, control, and state of the system to the output file
        output_file << this->m_t;
        for (double &value : u) {
            output_file << "," << value;
        }
        for (double &value : this->m_x) {
            output_file << "," << value;
        }
        output_file << std::endl;

        // If the time horizon stretches beyond the final time, then shrink it such that it doesn't
        if (this->m_t + horizon > time) {
            horizon = time - this->m_t;
        }

        // Update the control with a newly computed control input with the current state of the system
        int status = this->m_system_ptr->updateControl(this->m_x, this->m_x_star, this->m_t, horizon);
        if (status != 0) {
            return -1;
        }

        // Use odeint to integrate the system from the current time to the time of the next control update
        boost::numeric::odeint::integrate(*(this->m_system_ptr), this->m_x , this->m_t,
                                          this->m_t + this->m_mpc_time_step , this->m_mpc_time_step);

        // Update the time of the simulation
        this->m_t = this->m_t + this->m_mpc_time_step;
    }

    // Close the output file
    output_file.close();

    // Return with no errors
    return 0;
}