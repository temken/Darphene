#include "gtest/gtest.h"

#include <functional>
#include <limits>

#include "libphysica/Integration.hpp"
#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Natural_Units.hpp"

#include "graphene/Carbon_Wavefunctions.hpp"

using namespace libphysica::natural_units;
using namespace graphene;

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionValues)
{
	// ARRANGE
	Eigen::Vector3d position({Bohr_Radius, Bohr_Radius / 2, Bohr_Radius / 3});
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	double phi_2s  = -6.3755250780151725e-10;
	double phi_2px = 1.2651040375412864e-09;
	double phi_2py = 6.325520187706432e-10;
	double phi_2pz = 4.2170134584709545e-10;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position(position, "2s"), phi_2s);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position(position, "2px"), phi_2px);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position(position, "2py"), phi_2py);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position(position, "2pz"), phi_2pz);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization)
{
	// ARRANGE
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	std::vector<std::string> orbitals = {"2s", "2px", "2py", "2pz"};
	// ACT & ASSERT
	for(auto orbital : orbitals)
		EXPECT_FLOAT_EQ(hydrogenic.Normalization_Position(orbital), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestMomentumWavefunctionValues)
{
	// ARRANGE
	Eigen::Vector3d momentum({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3 / Bohr_Radius});
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	double phi_2s  = 1.5170029e7;
	double phi_2px = 2.33385e6;
	double phi_2py = 4.6677e6;
	double phi_2pz = 7.00155e6;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum(momentum, "2s")), phi_2s);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum(momentum, "2px")), phi_2px);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum(momentum, "2py")), phi_2py);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum(momentum, "2pz")), phi_2pz);
}

TEST(TestHydrogenicWavefunctions, TestMomentumWavefunctionNormalization)
{
	// ARRANGE
	double tol	= 1e-6;
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	std::vector<std::string> orbitals = {"2s", "2px", "2py", "2pz"};
	// ACT & ASSERT
	for(auto orbital : orbitals)
		EXPECT_NEAR(hydrogenic.Normalization_Momentum(orbital), 1.0, tol);
}