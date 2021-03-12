#ifndef __Multidimensional_Integration_hpp_
#define __Multidimensional_Integration_hpp_

#include <Eigen/Geometry>
#include <cmath>
#include <functional>
#include <string>

namespace graphene
{
extern Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi);

extern double Integrate_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, const std::string& method = "Trapezoidal");
extern double Integrate_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, const std::string& method = "Trapezoidal");
extern double Integrate_3D(std::function<double(Eigen::Vector3d)>& integrand, double rMin, double rMax, double cos_theta_min = -1.0, double cos_theta_max = 1.0, double phi_min = 0.0, double phi_max = 2.0 * M_PI, const std::string& method = "Trapezoidal");

}	// namespace graphene

#endif