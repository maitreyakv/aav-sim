/*
 * Copyright (C) 2020 Maitreya Venkataswamy - All Rights Reserved
 */

#include "SphereObstacle.h"

// TODO: Add documentation/comments

bool SphereObstacle::isLineCollide(arma::vec point_1, arma::vec point_2) {
    // https://en.wikipedia.org/wiki/Line-sphere_intersection

    // Obtain the unit vector defining the direction of the line
    arma::vec l = arma::normalise(point_2 - point_1);

    double term_1 = pow(arma::dot(l, (point_1 - this->m_center)), 2);
    double term_2 = pow(arma::norm(point_1 - this->m_center), 2) - pow(this->m_radius, 2);

    return term_1 - term_2 >= 0.0;
}