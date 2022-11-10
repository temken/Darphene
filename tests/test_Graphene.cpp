#include "gtest/gtest.h"

#include <cmath>
#include <random>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "graphene/Direct_Detection_Standard.hpp"
#include "graphene/Graphene.hpp"

using namespace graphene;
using namespace libphysica::natural_units;

TEST(TestGraphene, TestOverlapIntegralsHydrogenic)
{
	// ARRANGE
	Graphene graphene("Hydrogenic");
	std::vector<std::string> overlap_parameters = {"s", "sPrime", "Sss", "Ssigma", "Ssp"};
	std::vector<double> results					= {0.129, 0.00866, 0.212, 0.146, 0.159931};
	double tol									= 1.0e-3;
	// ACT & ASSERT
	for(unsigned int i = 0; i < overlap_parameters.size(); i++)
		EXPECT_NEAR(graphene.Overlap_Integral(overlap_parameters[i]), results[i], tol);
}

// TEST(TestGraphene, TestOverlapIntegralsRHF)
// {
// 	// ARRANGE
// 	Graphene graphene("RHF");
// 	std::vector<std::string> overlap_parameters = {"s", "sPrime", "Sss", "Ssigma", "Ssp"};
// 	std::vector<double> results					= {0.129, 0.00866, 0.212, 0.146, 0.159931};
// 	double tol									= 1.0e-3;
// 	// ACT & ASSERT
// 	for(unsigned int i = 0; i < overlap_parameters.size(); i++)
// 		EXPECT_NEAR(graphene.Overlap_Integral(overlap_parameters[i]), results[i], tol);
// }

TEST(TestGraphene, TestEnergyDispersionPi)
{
	// ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene("Hydrogenic");
	std::vector<std::vector<double>> results = {{-6.5602 * eV, 14.843393 * eV},
												{-2.6864482 * eV, 3.4822043 * eV},
												{0.0 * eV, 0.0 * eV}};
	double tolerance						 = 1.0e-5 * eV;
	// ACT & ASSERT
	int i = 0;
	for(auto& lVec : high_symmetry_points)
	{
		std::vector<double> energy_analytic = graphene.Energy_Dispersion_Pi_Analytic(lVec);
		EXPECT_NEAR(energy_analytic[0], results[i][0], tolerance);
		EXPECT_NEAR(energy_analytic[1], results[i++][1], tolerance);
	}
}

TEST(TestGraphene, TestEnergyDispersionPiAnalytic)
{
	// ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene("Hydrogenic");
	// ACT & ASSERT
	for(auto& lVec : high_symmetry_points)
	{
		std::vector<double> energy_analytic = graphene.Energy_Dispersion_Pi_Analytic(lVec);
		std::vector<double> energy_numeric	= graphene.Energy_Dispersion_Pi(lVec);
		EXPECT_FLOAT_EQ(energy_analytic[0], energy_numeric[0]);
		EXPECT_FLOAT_EQ(energy_analytic[1], energy_numeric[1]);
	}
}

TEST(TestGraphene, TestEnergyDispersionSigma)
{
	// ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene("Hydrogenic");
	std::vector<std::vector<double>> results = {
		{-17.833129584352082 * eV, -2.93125 * eV, -2.93125 * eV, 3.08466 * eV, 3.08466 * eV, 31.425824175824184 * eV},
		{-14.802126408956223 * eV, -11.691743504673768 * eV, -7.06817 * eV, 10.476314847203604 * eV, 12.661549197487779 * eV, 20.1867 * eV},
		{-13.064658876661905 * eV, -13.064658876661898 * eV, -8.56991 * eV, 12.119618979416404 * eV, 12.119618979416412 * eV, 20.604255319148937 * eV}};
	double tolerance = 1.0e-5 * eV;
	// ACT & ASSERT
	int i = 0;
	for(auto& lVec : high_symmetry_points)
	{
		std::vector<double> energy = graphene.Energy_Dispersion_Sigma(lVec);
		EXPECT_NEAR(energy[0], results[i][0], tolerance);
		EXPECT_NEAR(energy[1], results[i][1], tolerance);
		EXPECT_NEAR(energy[2], results[i][2], tolerance);
		EXPECT_NEAR(energy[3], results[i][3], tolerance);
		EXPECT_NEAR(energy[4], results[i][4], tolerance);
		EXPECT_NEAR(energy[5], results[i++][5], tolerance);
	}
}

TEST(TestGraphene, TestEnergyBands)
{
	// ARRANGE
	Graphene graphene("Hydrogenic");
	int N			 = 100;
	double tolerance = 1e-22;

	// ACT
	std::vector<std::vector<double>> energy_bands = graphene.Energy_Bands(N);

	// ASSERT
	EXPECT_EQ(energy_bands.size(), N);
	for(auto& entry : energy_bands)
		EXPECT_EQ(entry.size(), 9);
	for(int i = 1; i < 9; i++)
		EXPECT_NEAR(energy_bands[0][i], energy_bands[N - 1][i], tolerance);
}

TEST(TestGraphene, TestBZ)
{
	// ARRANGE
	Graphene graphene("Hydrogenic");
	std::random_device rd;
	std::mt19937 PRNG(rd());
	double kMax = 20. * graphene.b;
	double tol	= 1e-10;
	//
	// ACT & ASSERT
	for(int i = 0; i < 100; i++)
	{
		Eigen::Vector3d k = {libphysica::Sample_Uniform(PRNG, -kMax, kMax), libphysica::Sample_Uniform(PRNG, -kMax, kMax), 0.0};
		Eigen::Vector3d G = graphene.Find_G_Vector(k);
		Eigen::Vector3d l = k - G;
		EXPECT_TRUE(graphene.In_1BZ(l));
		// Check that coefficients are integer.
		double m = 1.0 / graphene.b * (G[0] + G[1] / sqrt(3.0));
		double n = 1.0 / graphene.b * (G[0] - G[1] / sqrt(3.0));
		EXPECT_NEAR(m, std::round(m), tol);
		EXPECT_NEAR(n, std::round(n), tol);
	}
}

TEST(TestGraphene, TestResponseFunctionNormalizationHydrogenic)
{
	// ARRANGE
	double tol = 0.1;
	Graphene graphene("Hydrogenic");
	std::function<double(double, double, double)> integrand = [&graphene](double l, double cos_theta, double phi) {
		Eigen::Vector3d lVec = Spherical_Coordinates(l, acos(cos_theta), phi);
		double W			 = 0.0;
		for(int band = 0; band < 4; band++)
			W += graphene.Material_Response_Function(band, lVec);
		return l * l * W;
	};
	double kMax = 25.0 * keV;
	// ACT
	double norm = libphysica::Integrate_3D(integrand, 0, kMax, -1.0, 1.0, 0.0, 2 * M_PI, "Vegas", 1000);
	// ASSERT
	ASSERT_NEAR(norm, 4.0, tol);
}

TEST(TestGraphene, TestResponseFunctionNormalizationRHF)
{
	// ARRANGE
	double tol = 0.1;
	Graphene graphene("RHF");
	std::function<double(double, double, double)> integrand = [&graphene](double l, double cos_theta, double phi) {
		Eigen::Vector3d lVec = Spherical_Coordinates(l, acos(cos_theta), phi);
		double W			 = 0.0;
		for(int band = 0; band < 4; band++)
			W += graphene.Material_Response_Function(band, lVec);
		return l * l * W;
	};
	double kMax = 25.0 * keV;
	// ACT
	double norm = libphysica::Integrate_3D(integrand, 0, kMax, -1.0, 1.0, 0.0, 2 * M_PI, "Vegas", 1000);
	// ASSERT
	ASSERT_NEAR(norm, 4.0, tol);
}