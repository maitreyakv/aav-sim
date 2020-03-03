/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

#include "QuadcopterDynamics.h"
#include "../../mpcsim/src/DDPController.h"
#include "../../mpcsim/src/QuadraticCost.h"
#include "../../mpcsim/src/Simulation.h"

arma::vec genVectorFromInputOption(boost::property_tree::ptree inputs, std::string option) {
    std::string str = inputs.get<std::string>(option);
    std::vector<std::string> str_split;
    boost::split(str_split , str, boost::is_any_of(","));

    arma::vec vector = arma::zeros<arma::vec>( str_split.size() );

    for (long unsigned int i = 0; i < str_split.size(); i++) {
        vector(i) = stod(str_split[i]);
    }

    return vector;
}

arma::mat genDiagonalMatrixFromInputOption(boost::property_tree::ptree inputs, std::string option) {
    std::string str = inputs.get<std::string>(option);
    std::vector<std::string> str_split;
    boost::split(str_split , str, boost::is_any_of(","));

    arma::mat matrix = arma::zeros<arma::mat>( str_split.size(), str_split.size() );

    for (long unsigned int i = 0; i < str_split.size(); i++) {
        matrix(i,i) = stod(str_split[i]);
    }

    return matrix;
}

int main() {

    boost::property_tree::ptree inputs;
    boost::property_tree::ini_parser::read_ini("example_input.ini", inputs);

    double m    = stod( inputs.get<std::string>("QuadcopterParameters.m")    );
    double g    = stod( inputs.get<std::string>("QuadcopterParameters.g")    );
    double l    = stod( inputs.get<std::string>("QuadcopterParameters.l")    );
    double k    = stod( inputs.get<std::string>("QuadcopterParameters.k")    );
    double b    = stod( inputs.get<std::string>("QuadcopterParameters.b")    );
    double I_xx = stod( inputs.get<std::string>("QuadcopterParameters.I_xx") );
    double I_yy = stod( inputs.get<std::string>("QuadcopterParameters.I_yy") );
    double I_zz = stod( inputs.get<std::string>("QuadcopterParameters.I_zz") );

    QuadcopterDynamics* quadcopter_dynamics_ptr = new QuadcopterDynamics(m, g, l, k, b, I_xx, I_yy, I_zz);

    arma::vec x_0 = genVectorFromInputOption(inputs, "SimulationParameters.x_0");

    arma::vec x_star = genVectorFromInputOption(inputs, "SimulationParameters.x_star");

    arma::mat Q_f = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.Q_f");
    arma::mat R = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.R");
    QuadraticCost* quadratic_cost_ptr = new QuadraticCost(Q_f, R);

    arma::vec u_max = stod( inputs.get<std::string>("QuadcopterParameters.u_max") ) * arma::ones<arma::vec>(4);

    DDPController* ddp_controller_ptr = new DDPController( quadcopter_dynamics_ptr, quadratic_cost_ptr, u_max,
        stod( inputs.get<std::string>("ControllerParameters.num_iterations") ),
        stod( inputs.get<std::string>("ControllerParameters.num_discretizations") ),
        stod( inputs.get<std::string>("ControllerParameters.learning_rate") ) );

    System* system_ptr = new System(quadcopter_dynamics_ptr, ddp_controller_ptr);

    Simulation* simulation_ptr = new Simulation(system_ptr,
        1.0 / stod( inputs.get<std::string>("SimulationParameters.mpc_rate") ),
        x_0, x_star, stod( inputs.get<std::string>("SimulationParameters.horizon") ));

    simulation_ptr->simulate( stod( inputs.get<std::string>("SimulationParameters.sim_time") ) );

    return 0;
}