#include "gtest/gtest.h"

#include <cmath>

#include "Multidimensional_Integration.hpp"

using namespace graphene;

TEST(TestMultidimensionalIntegration, TestTrapezoidal2D)
{
	// ARRANGE
	std::function<double(double, double)> f = [](double x, double y) {
		return x * y * y;
	};
	double result	   = 1.0 / 6.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	std::string method = "Trapezoidal";
	// ACT
	double integral = Integrate_2D(f, xMin, xMax, yMin, yMax, method);
	// ASSERT
	ASSERT_FLOAT_EQ(integral, result);
}

TEST(TestMultidimensionalIntegration, TestGaussKronrod2D)
{
	// ARRANGE
	std::function<double(double, double)> f = [](double x, double y) {
		return x * y * y;
	};
	double result	   = 1.0 / 6.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	std::string method = "Gauss-Kronrod";
	// ACT
	double integral = Integrate_2D(f, xMin, xMax, yMin, yMax, method);
	// ASSERT
	ASSERT_FLOAT_EQ(integral, result);
}

TEST(TestMultidimensionalIntegration, TestGauss2D)
{
	// ARRANGE
	std::function<double(double, double)> f = [](double x, double y) {
		return x * y * y;
	};
	double result	   = 1.0 / 6.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	std::string method = "Gauss";
	// ACT
	double integral = Integrate_2D(f, xMin, xMax, yMin, yMax, method);
	// ASSERT
	ASSERT_FLOAT_EQ(integral, result);
}

TEST(TestMultidimensionalIntegration, TestTrapezoidal3D)
{
	// ARRANGE
	std::function<double(double, double, double)> f = [](double x, double y, double z) {
		return x * y * pow(z, 6);
	};
	double result	   = 32.0 / 7.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	double zMin		   = 0.0;
	double zMax		   = 2.0;
	std::string method = "Trapezoidal";
	double tolerance   = 1e-5;
	// ACT
	double integral = Integrate_3D(f, xMin, xMax, yMin, yMax, zMin, zMax, method);
	// ASSERT
	ASSERT_NEAR(integral, result, tolerance);
}

TEST(TestMultidimensionalIntegration, TestGauss3D)
{
	// ARRANGE
	std::function<double(double, double, double)> f = [](double x, double y, double z) {
		return x * y * pow(z, 6);
	};
	double result	   = 32.0 / 7.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	double zMin		   = 0.0;
	double zMax		   = 2.0;
	std::string method = "Gauss";
	double tolerance   = 1e-5;
	// ACT
	double integral = Integrate_3D(f, xMin, xMax, yMin, yMax, zMin, zMax, method);
	// ASSERT
	ASSERT_NEAR(integral, result, tolerance);
}

TEST(TestMultidimensionalIntegration, TestGaussKronrod3D)
{
	// ARRANGE
	std::function<double(double, double, double)> f = [](double x, double y, double z) {
		return x * y * pow(z, 6);
	};
	double result	   = 32.0 / 7.0;
	double xMin		   = 0.0;
	double xMax		   = 1.0;
	double yMin		   = 0.0;
	double yMax		   = 1.0;
	double zMin		   = 0.0;
	double zMax		   = 2.0;
	std::string method = "Gauss-Kronrod";
	double tolerance   = 1e-5;
	// ACT
	double integral = Integrate_3D(f, xMin, xMax, yMin, yMax, zMin, zMax, method);
	// ASSERT
	ASSERT_NEAR(integral, result, tolerance);
}