/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <boost/numeric/odeint.hpp>

#include "Simulation.h"

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

void Simulation::simulate(double time) {
    // Make armadillo vectors resizable in odeint
    BOOST_STATIC_ASSERT(boost::numeric::odeint::is_resizeable<arma::vec>::value == true);

    // Initialize the time horizon with the specified horizon, this may be edited during the simulation
    double horizon = this->m_horizon;

    // Perform the simulation by updating the control and integrating the system
    while (this->m_t < time) {
        // TEMP: Temporary output for post-processing, will be replaced with proper file IO
        std::cout << this->m_t << "," << this->m_system_ptr->getControl()[0] << "," << this->m_x[0] << "," << this->m_x[2] << std::endl;

        // If the time horizon stretches beyond the final time, then shrink it such that it doesn't
        if (this->m_t + horizon > time) {
            horizon = time - this->m_t;
        }

        // Update the control with a newly computed control input with the current state of the system
        this->m_system_ptr->updateControl(this->m_x, this->m_x_star, this->m_t, horizon);

        // Use odeint to integrate the system from the current time to the time of the next control update
        boost::numeric::odeint::integrate(*(this->m_system_ptr), this->m_x , this->m_t ,
                                          this->m_t + this->m_mpc_time_step , this->m_mpc_time_step);

        // Update the time of the simulation
        this->m_t = this->m_t + this->m_mpc_time_step;
    }
}