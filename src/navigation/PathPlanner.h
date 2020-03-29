/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAV_SIM_PATHPLANNER_H
#define AAV_SIM_PATHPLANNER_H

// TODO: Add comments and documentation
#include <armadillo>

#include "../simulation/Obstacle.h"

class PathPlanner {

private:
    // List of obstacles to avoid
    std::vector<Obstacle*> m_obstacle_ptrs;

public:
    PathPlanner(std::vector<Obstacle*> obstacle_ptrs)
        : m_obstacle_ptrs(obstacle_ptrs)
    {}

    virtual std::vector<arma::vec> computePath(arma::vec pos_start, arma::vec pos_target) = 0;
};


#endif //AAV_SIM_PATHPLANNER_H
