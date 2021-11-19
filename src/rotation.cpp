#include "rotation.h"

std::vector<double> RotationX(double theta, std::vector<double> vec) {
    return {vec[0], vec[1]*cos(theta)-vec[2]*sin(theta), vec[1]*sin(theta)+vec[2]*cos(theta)};
}

std::vector<double> RotationY(double theta, std::vector<double> vec) {
    return {vec[0]*cos(theta)+vec[2]*sin(theta), vec[1], vec[2]*cos(theta)-vec[0]*sin(theta)};
}

std::vector<double> RotationZ(double theta, std::vector<double> vec) {
    return {vec[0]*cos(theta)-vec[1]*sin(theta), vec[0]*sin(theta)+vec[1]*cos(theta), vec[2]};
}

std::vector<double> Rotation(std::vector<double> angle, std::vector<double> vec) {
    return RotationZ(angle[2], RotationY(angle[1], RotationX(angle[0], vec)));
}