#include "gtest/gtest.h"

// Headers from libphysica
#include "Natural_Units.hpp"

#include "Hydrogenic_Wavefunctions.hpp"

using namespace libphysica::natural_units;
using namespace graphene;

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionValues)
{
	//ARRANGE
	Eigen::Vector3d position({Bohr_Radius, Bohr_Radius / 2, Bohr_Radius / 3});
	double Zeff	   = 2.0;
	double phi_2s  = -6.3755250780151725e-10;
	double phi_2px = 1.2651040375412864e-09;
	double phi_2py = 6.325520187706432e-10;
	double phi_2pz = 4.2170134584709545e-10;
	//ACT & ASSERT
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_2s(position, Zeff), phi_2s);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_2px(position, Zeff), phi_2px);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_2py(position, Zeff), phi_2py);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_2pz(position, Zeff), phi_2pz);
}

// TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization)
// {
// 	//ARRANGE
// 	//ACT & ASSERT

// }

TEST(TestHydrogenicWavefunctions, TestMomentumWavefunctionValues)
{
	//ARRANGE
	Eigen::Vector3d momentum({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3 / Bohr_Radius});
	double Zeff	   = 2.0;
	double phi_2s  = 1.5170029e7;
	double phi_2px = 2.33385e6;
	double phi_2py = 4.6677e6;
	double phi_2pz = 7.00155e6;
	//ACT & ASSERT
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum_2s(momentum, Zeff), phi_2s);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum_2px(momentum, Zeff), phi_2px);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum_2py(momentum, Zeff), phi_2py);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum_2pz(momentum, Zeff), phi_2pz);
}