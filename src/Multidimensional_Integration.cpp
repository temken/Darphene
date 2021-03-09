#include "Multidimensional_Integration.hpp"

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <iostream>

namespace graphene
{
using namespace boost::math::quadrature;

// 1. Two-dimensional integrals
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

// 2. Three-dimensional integrals
double Trapezoidal_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
	auto integrand_x = [&integrand, yMin, yMax, zMin, zMax](double x) {
		auto integrand_y = [&integrand, x, zMin, zMax](double y) {
			auto integrand_z = [&integrand, x, y](double z) {
				return integrand(x, y, z);
			};
			double integral_z = trapezoidal(integrand_z, zMin, zMax);
			return integral_z;
		};
		double integral_y = trapezoidal(integrand_y, yMin, yMax);
		return integral_y;
	};
	return trapezoidal(integrand_x, xMin, xMax);
}

double Gauss_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
	auto integrand_x = [&integrand, yMin, yMax, zMin, zMax](double x) {
		auto integrand_y = [&integrand, x, zMin, zMax](double y) {
			auto integrand_z = [&integrand, x, y](double z) {
				return integrand(x, y, z);
			};
			double integral_z = gauss<double, 30>::integrate(integrand_z, zMin, zMax);
			return integral_z;
		};
		double integral_y = gauss<double, 30>::integrate(integrand_y, yMin, yMax);
		return integral_y;
	};
	return gauss<double, 30>::integrate(integrand_x, xMin, xMax);
}

double Gauss_Kronrod_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
	auto integrand_x = [&integrand, yMin, yMax, zMin, zMax](double x) {
		auto integrand_y = [&integrand, x, zMin, zMax](double y) {
			auto integrand_z = [&integrand, x, y](double z) {
				return integrand(x, y, z);
			};
			double integral_z = gauss_kronrod<double, 31>::integrate(integrand_z, zMin, zMax, 5, 1e-9);
			return integral_z;
		};
		double integral_y = gauss_kronrod<double, 31>::integrate(integrand_y, yMin, yMax, 5, 1e-9);
		return integral_y;
	};
	return gauss_kronrod<double, 31>::integrate(integrand_x, xMin, xMax, 5, 1e-9);
}

double Integrate_3D(std::function<double(double, double, double)>& integrand, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax, const std::string& method)
{
	if(method == "Trapezoidal")
		return Trapezoidal_3D(integrand, xMin, xMax, yMin, yMax, zMin, zMax);
	else if(method == "Gauss")
		return Gauss_3D(integrand, xMin, xMax, yMin, yMax, zMin, zMax);
	else if(method == "Gauss-Kronrod")
		return Gauss_Kronrod_3D(integrand, xMin, xMax, yMin, yMax, zMin, zMax);
	else
	{
		std::cerr << "Error in Integrate_3D(): Method " << method << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace graphene
