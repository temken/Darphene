#include "Multidimensional_Integration.hpp"

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <iostream>

namespace graphene
{
using namespace boost::math::quadrature;

double Trapezoidal_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax)
{
	auto integrand_x = [&integrand, yMin, yMax](double x) {
		auto integrand_y = [&integrand, x](double y) {
			return integrand(x, y);
		};
		double integral_y = trapezoidal(integrand_y, yMin, yMax);
		return integral_y;
	};
	return trapezoidal(integrand_x, xMin, xMax);
}

double Gauss_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax)
{
	auto integrand_x = [&integrand, yMin, yMax](double x) {
		auto integrand_y = [&integrand, x](double y) {
			return integrand(x, y);
		};
		double integral_y = gauss<double, 30>::integrate(integrand_y, yMin, yMax);
		return integral_y;
	};
	return gauss<double, 30>::integrate(integrand_x, xMin, xMax);
}

double Gauss_Kronrod_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax)
{
	auto integrand_x = [&integrand, yMin, yMax](double x) {
		auto integrand_y = [&integrand, x](double y) {
			return integrand(x, y);
		};
		double integral_y = gauss_kronrod<double, 31>::integrate(integrand_y, yMin, yMax, 5, 1e-9);
		return integral_y;
	};
	return gauss_kronrod<double, 31>::integrate(integrand_x, xMin, xMax, 5, 1e-9);
}

double Integrate_2D(std::function<double(double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, const std::string& method)
{
	if(method == "Trapezoidal")
		return Trapezoidal_2D(integrand, xMin, xMax, yMin, yMax);
	else if(method == "Gauss")
		return Gauss_2D(integrand, xMin, xMax, yMin, yMax);
	else if(method == "Gauss-Kronrod")
		return Gauss_Kronrod_2D(integrand, xMin, xMax, yMin, yMax);
	else
	{
		std::cerr << "Error in Integrate_2D(): Method " << method << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace graphene
