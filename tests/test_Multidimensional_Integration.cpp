#include "gtest/gtest.h"

#include "Multidimensional_Integration.hpp"

using namespace graphene;

TEST(TestMultidimensionalIntegration, TestTrapezoidal)
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

TEST(TestMultidimensionalIntegration, TestGaussKronrod)
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

TEST(TestMultidimensionalIntegration, TestGauss)
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