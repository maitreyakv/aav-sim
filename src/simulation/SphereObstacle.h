/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#ifndef AAV_SIM_SPHEREOBSTACLE_H
#define AAV_SIM_SPHEREOBSTACLE_H

// TODO: add documentation/comments
#include "Obstacle.h"

class SphereObstacle : public Obstacle {
private:
    // Center of the sphere
    arma::vec m_center;

    // Radius of the sphere
    double m_radius;
    
    // Factor to multiply the actual radius to get the perceived radius
    double m_increase_factor;
    
public:
    SphereObstacle(arma::vec center, double radius, double increase_factor)
        : m_center(center),
          m_radius(radius),
          m_increase_factor(increase_factor)
    {}

    // https://en.wikipedia.org/wiki/Line-sphere_intersection
    bool isLineCollide(arma::vec point_1, arma::vec point_2);
};


#endif //AAV_SIM_SPHEREOBSTACLE_H
