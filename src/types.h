//
// Created by Gerrit Freiwald on 08.05.2023.
//

#include <Eigen/Dense>

#ifndef MY_MD_CODE_TYPES_H
#define MY_MD_CODE_TYPES_H

typedef Eigen::Array3Xd Positions_t;
typedef Eigen::Array3Xd Velocities_t;
typedef Eigen::Array3Xd Forces_t;
typedef Eigen::ArrayXd  Masses_t;
typedef Eigen::Array3Xd Acceleration_t;
typedef Eigen::ArrayXd  Potential_t;
typedef std::vector<std::string> Names_t;


#endif // MY_MD_CODE_TYPES_H
