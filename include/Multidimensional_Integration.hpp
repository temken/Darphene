#ifndef __Multidimensional_Integration_hpp_
#define __Multidimensional_Integration_hpp_

#include <functional>
#include <string>

namespace graphene
{

extern double Integrate_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, const std::string& method = "Trapezoidal");
extern double Integrate_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, const std::string& method = "Trapezoidal");

}	// namespace graphene

#endif