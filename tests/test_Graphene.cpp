#include "gtest/gtest.h"

// Headers from libphysica
#include "Natural_Units.hpp"

#include "Graphene.hpp"

using namespace graphene;
using namespace libphysica::natural_units;

TEST(TestGraphene, TestEnergyDispersionPi)
{
	//ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene;
	std::vector<std::vector<double>> results = {{-6.5602 * eV, 14.843393 * eV},
												{-2.6864482 * eV, 3.4822043 * eV},
												{0.0 * eV, 0.0 * eV}};
	double tolerance						 = 1.0e-5 * eV;
	//ACT & ASSERT
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
	//ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene;
	//ACT & ASSERT
	for(auto& lVec : high_symmetry_points)
	{
		std::vector<double> energy_analytic = graphene.Energy_Dispersion_Pi_Analytic(lVec);
		std::vector<double> energy_numeric	= graphene.Energy_Dispersion_Pi(lVec);
		EXPECT_FLOAT_EQ(energy_analytic[0], energy_numeric[0]);
		EXPECT_FLOAT_EQ(energy_analytic[1], energy_numeric[1]);
	}
}