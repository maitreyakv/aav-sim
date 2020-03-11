#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

#include "CartPoleDynamics.h"
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

/**
 * Main point of entry for program, sets up the simulation and executes the simulation
 *
 * @return 0 if successful, nonzero otherwise
 */
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

    // Get the system parameters from the input file
    double mass_cart  = stod( inputs.get<std::string>("CartPoleParameters.mass_cart") );
    double mass_pole  = stod( inputs.get<std::string>("CartPoleParameters.mass_pole") );
    double length_pole  = stod( inputs.get<std::string>("CartPoleParameters.length_pole") );

    // Get errors of system parameters from input file
    double error_mass_cart  = stod( inputs.get<std::string>("CartPoleParameters.error_mass_cart") );
    double error_mass_pole  = stod( inputs.get<std::string>("CartPoleParameters.error_mass_pole") );
    double error_length_pole  = stod( inputs.get<std::string>("CartPoleParameters.error_length_pole") );

    // Initialize the true and false Dynamics of the System
    CartPoleDynamics* cart_pole_dynamics_true_ptr = new CartPoleDynamics(mass_cart, mass_pole, length_pole);
    CartPoleDynamics* cart_pole_dynamics_false_ptr = new CartPoleDynamics(mass_cart + error_mass_cart,
                                                                          mass_pole + error_mass_pole,
                                                                          length_pole + error_length_pole);

    // Read the cost matrices from the input file
    arma::mat Q_f = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.Q_f");
    arma::mat R = genDiagonalMatrixFromInputOption(inputs, "CostFunctionParameters.R");

    // Initialize the Quadratic cost functions
    QuadraticCost* quadratic_cost_ptr = new QuadraticCost(Q_f, R);

    // Read the maximum control from the input file
    arma::vec u_max = genVectorFromInputOption(inputs, "CartPoleParameters.u_max");

    // Initialize the DDPController
    DDPController* ddp_controller_ptr = new DDPController(cart_pole_dynamics_false_ptr, quadratic_cost_ptr, u_max,
        stod( inputs.get<std::string>("ControllerParameters.num_discretizations") ),
        stod( inputs.get<std::string>("ControllerParameters.num_iterations") ),
        stod( inputs.get<std::string>("ControllerParameters.learning_rate") ) );

    // Initialize the System
    System* system_ptr = new System(cart_pole_dynamics_true_ptr, ddp_controller_ptr);

    // Read the initial state from the input file
    arma::vec x_0 = genVectorFromInputOption(inputs, "SimulationParameters.x_0");

    // Read the target state from the input file
    arma::vec x_star = genVectorFromInputOption(inputs, "SimulationParameters.x_star");

    // Initialize the Simulation
    Simulation* simulation_ptr = new Simulation(system_ptr,
        1.0 / stod( inputs.get<std::string>("SimulationParameters.mpc_rate") ),
        x_0, x_star, stod( inputs.get<std::string>("SimulationParameters.horizon") ));

    // Execute the simulation
    int status = simulation_ptr->simulate( stod( inputs.get<std::string>("SimulationParameters.sim_time") ) );

    // Return with no errors
    return status;
}