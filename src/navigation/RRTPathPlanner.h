/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAV_SIM_RRTPATHPLANNER_H
#define AAV_SIM_RRTPATHPLANNER_H

// TODO: add documentation/comments
#include "PathPlanner.h"

class RRTPathPlanner : public PathPlanner {
public:
    RRTPathPlanner(std::vector<Obstacle*> obstacle_ptrs)
        : PathPlanner(obstacle_ptrs)
    {}

    std::vector<arma::vec> computePath(arma::vec pos_start, arma::vec pos_target);

};


#endif //AAV_SIM_RRTPATHPLANNER_H
