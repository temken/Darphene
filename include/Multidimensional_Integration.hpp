#ifndef __Multidimensional_Integration_hpp_
#define __Multidimensional_Integration_hpp_

#include <functional>
#include <string>

// #include <boost/math/quadrature/exp_sinh.hpp>
// #include <boost/math/quadrature/gauss.hpp>
// #include <boost/math/quadrature/gauss_kronrod.hpp>
// #include <boost/math/quadrature/naive_monte_carlo.hpp>
// #include <boost/math/quadrature/trapezoidal.hpp>

namespace graphene
{

extern double Integrate_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, const std::string& method = "Trapezoidal");

}	// namespace graphene

#endif