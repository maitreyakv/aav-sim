/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "SphereObstacle.h"

// TODO: Add documentation/comments

bool SphereObstacle::isLineCollide(arma::vec point_1, arma::vec point_2) {
    // https://cseweb.ucsd.edu/classes/sp19/cse291-d/Files/CSE291_13_CollisionDetection.pdf

    double t = arma::dot( (this->m_center - point_1), (point_2 - point_1) / pow(arma::norm(point_2 - point_1), 2) );

    bool is_collide;

    if (t <= 0.0) {
        is_collide = pow(arma::norm(this->m_center - point_1), 2) <= pow(this->m_increase_factor * this->m_radius, 2);
    } else if (t >= 1.0) {
        is_collide = pow(arma::norm(this->m_center - point_2), 2) <= pow(this->m_increase_factor * this->m_radius, 2);
    } else {
        arma::vec x = point_1 + t * (point_2 - point_1);
        is_collide = pow(arma::norm(this->m_center - x), 2) <= pow(this->m_increase_factor * this->m_radius, 2);
    }

    return is_collide;
}