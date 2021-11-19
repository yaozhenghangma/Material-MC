#ifndef ROTATION
#define ROTATION

#include <vector>
#include <math.h>

std::vector<double> RotationX(double theta, std::vector<double> vec);
std::vector<double> RotationY(double theta, std::vector<double> vec);
std::vector<double> RotationZ(double theta, std::vector<double> vec);
std::vector<double> Rotation(std::vector<double> angle, std::vector<double> vec);

#endif