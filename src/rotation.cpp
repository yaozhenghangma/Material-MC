#include "rotation.h"

/**
 * @file rotation.cpp
 * @brief 3D vector rotation helpers used by spin initialization and updates.
 */

/**
 * @brief Rotates a 3D vector around the X axis.
 *
 * @param theta Rotation angle in radians.
 * @param vec Input vector in Cartesian coordinates: {x, y, z}.
 * @return std::vector<double> Rotated vector {x', y', z'}.
 */
std::vector<double> RotationX(double theta, std::vector<double> vec) {
    return {vec[0], vec[1]*cos(theta)-vec[2]*sin(theta), vec[1]*sin(theta)+vec[2]*cos(theta)};
}

/**
 * @brief Rotates a 3D vector around the Y axis.
 *
 * @param theta Rotation angle in radians.
 * @param vec Input vector in Cartesian coordinates: {x, y, z}.
 * @return std::vector<double> Rotated vector {x', y', z'}.
 */
std::vector<double> RotationY(double theta, std::vector<double> vec) {
    return {vec[0]*cos(theta)+vec[2]*sin(theta), vec[1], vec[2]*cos(theta)-vec[0]*sin(theta)};
}

/**
 * @brief Rotates a 3D vector around the Z axis.
 *
 * @param theta Rotation angle in radians.
 * @param vec Input vector in Cartesian coordinates: {x, y, z}.
 * @return std::vector<double> Rotated vector {x', y', z'}.
 */
std::vector<double> RotationZ(double theta, std::vector<double> vec) {
    return {vec[0]*cos(theta)-vec[1]*sin(theta), vec[0]*sin(theta)+vec[1]*cos(theta), vec[2]};
}

/**
 * @brief Applies Euler-style axis rotations to a vector in X->Y->Z order.
 *
 * The function first rotates around X by angle[0], then around Y by angle[1],
 * and finally around Z by angle[2].
 *
 * @param angle Rotation angles in radians: {theta_x, theta_y, theta_z}.
 * @param vec Input vector in Cartesian coordinates: {x, y, z}.
 * @return std::vector<double> Rotated vector after sequential X, Y, Z rotations.
 */
std::vector<double> Rotation(std::vector<double> angle, std::vector<double> vec) {
    return RotationZ(angle[2], RotationY(angle[1], RotationX(angle[0], vec)));
}