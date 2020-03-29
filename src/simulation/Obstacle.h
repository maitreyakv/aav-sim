/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAV_SIM_OBSTACLE_H
#define AAV_SIM_OBSTACLE_H

// TODO: Add documentation/comments
#include <armadillo>

class Obstacle {
public:
    virtual bool isLineCollide(arma::vec point_1, arma::vec point1) = 0;
};


#endif //AAV_SIM_OBSTACLE_H
