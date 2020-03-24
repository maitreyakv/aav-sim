/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

#include "QuadcopterDynamics.h"
#include "../../mpcsim/src/DDPController.h"
#include "../../mpcsim/src/SamplerController.h"
#include "../../mpcsim/src/QuadraticCost.h"
#include "../../mpcsim/src/QuadraticObstacleCost.h"
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

int main(int argc, char** argv) {
    // Check if all program inputs were given
    if (argc < 2) {
        std::cout << "error: not enough program arguments provided\nrun program as:\n";
        std::cout << "    $ " << argv[0] << " [input_file]" << std::endl;
        return -1;
    }

    // Get the name of the configuration file
    std::string input_file_name = argv[1];

    // Read the configuration file with the simulation inputs
    boost::property_tree::ptree inputs;

    // Attempt to read the input file
    try {
        boost::property_tree::ini_parser::read_ini(input_file_name, inputs);
    } catch (boost::wrapexcept<boost::property_tree::ini_parser::ini_parser_error> &e) {
        std::cout << "error: problem reading input file\n";
        return -1;
    }

    double m    = stod( inputs.get<std::string>("QuadcopterParameters.m")    );
    double g    = stod( inputs.get<std::string>("QuadcopterParameters.g")    );
    double l    = stod( inputs.get<std::string>("QuadcopterParameters.l")    );
    double k    = stod( inputs.get<std::string>("QuadcopterParameters.k")    );
    double b    = stod( inputs.get<std::string>("QuadcopterParameters.b")    );
    double I_xx = stod( inputs.get<std::string>("QuadcopterParameters.I_xx") );
    double I_yy = stod( inputs.get<std::string>("QuadcopterParameters.I_yy") );
    double I_zz = stod( inputs.get<std::string>("QuadcopterParameters.I_zz") );

    QuadcopterDynamics* dynamics_ptr = new QuadcopterDynamics(m, g, l, k, b, I_xx, I_yy, I_zz);

    arma::vec x_0 = genVectorFromInputOption(inputs, "SimulationParameters.x_0");

    arma::vec x_star = genVectorFromInputOption(inputs, "SimulationParameters.x_star");

    arma::mat Q_f = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.Q_f");
    arma::mat R = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.R");
    std::vector<arma::vec> obstacles;
    obstacles.push_back(arma::zeros<arma::vec>(12));
    arma::mat sigma = arma::zeros<arma::mat>(12, 12);
    sigma(0,0) = stod( inputs.get<std::string>("CostFunctionParameters.sigma") );
    sigma(1,1) = sigma(0,0);
    sigma(2,2) = sigma(0,0);
    QuadraticObstacleCost* cost_ptr = new QuadraticObstacleCost(Q_f, R, obstacles, sigma);

    arma::vec u_max = stod( inputs.get<std::string>("QuadcopterParameters.u_max") ) * arma::ones<arma::vec>(4);

    DDPController* controller_ptr = new DDPController( dynamics_ptr, cost_ptr, u_max,
        stod( inputs.get<std::string>("DDPControllerParameters.num_discretizations") ),
        stod( inputs.get<std::string>("DDPControllerParameters.num_iterations") ),
        stod( inputs.get<std::string>("DDPControllerParameters.learning_rate") ) );

    //SamplerController* controller_ptr = new SamplerController( dynamics_ptr, cost_ptr, u_max,
    //    stod( inputs.get<std::string>("SamplerControllerParameters.num_discretizations") ),
    //    stod( inputs.get<std::string>("SamplerControllerParameters.num_samples") ) );

    System* system_ptr = new System(dynamics_ptr, controller_ptr);

    Simulation* simulation_ptr = new Simulation(system_ptr,
        1.0 / stod( inputs.get<std::string>("SimulationParameters.mpc_rate") ),
        x_0, x_star, stod( inputs.get<std::string>("SimulationParameters.horizon") ));

    simulation_ptr->simulate( stod( inputs.get<std::string>("SimulationParameters.sim_time") ) );

    return 0;
}