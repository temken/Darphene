#include "gtest/gtest.h"

#include <functional>
#include <limits>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Numerics.hpp"

#include "Hydrogenic_Wavefunctions.hpp"
#include "Multidimensional_Integration.hpp"

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
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction(position, "2s", Zeff), phi_2s);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction(position, "2px", Zeff), phi_2px);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction(position, "2py", Zeff), phi_2py);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction(position, "2pz", Zeff), phi_2pz);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2s)
{
	//ARRANGE
	double Zeff										 = 2.0;
	std::function<double(Eigen::Vector3d)> integrand = [Zeff](Eigen::Vector3d rVec) {
		double phi = Hydrogenic_Wavefunction(rVec, "2s", Zeff);
		return phi * phi;
	};

	//ACT & ASSERT
	ASSERT_FLOAT_EQ(Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI, "Gauss"), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2px)
{
	//ARRANGE
	double Zeff										 = 2.0;
	std::function<double(Eigen::Vector3d)> integrand = [Zeff](Eigen::Vector3d rVec) {
		double phi = Hydrogenic_Wavefunction(rVec, "2px", Zeff);
		return phi * phi;
	};

	//ACT & ASSERT
	ASSERT_FLOAT_EQ(Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI, "Gauss"), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2py)
{
	//ARRANGE
	double Zeff										 = 2.0;
	std::function<double(Eigen::Vector3d)> integrand = [Zeff](Eigen::Vector3d rVec) {
		double phi = Hydrogenic_Wavefunction(rVec, "2py", Zeff);
		return phi * phi;
	};

	//ACT & ASSERT
	ASSERT_FLOAT_EQ(Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI, "Gauss"), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2pz)
{
	//ARRANGE
	double Zeff										 = 2.0;
	std::function<double(Eigen::Vector3d)> integrand = [Zeff](Eigen::Vector3d rVec) {
		double phi = Hydrogenic_Wavefunction(rVec, "2pz", Zeff);
		return phi * phi;
	};

	//ACT & ASSERT
	ASSERT_FLOAT_EQ(Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI, "Gauss"), 1.0);
}

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
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum(momentum, "2s", Zeff), phi_2s);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum(momentum, "2px", Zeff), phi_2px);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum(momentum, "2py", Zeff), phi_2py);
	EXPECT_FLOAT_EQ(Hydrogenic_Wavefunction_Momentum(momentum, "2pz", Zeff), phi_2pz);
}