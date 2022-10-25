#include "gtest/gtest.h"

#include <functional>
#include <limits>

#include "libphysica/Integration.hpp"
#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Natural_Units.hpp"

#include "graphene/Hydrogenic_Wavefunctions.hpp"

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
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position_2s(position), phi_2s);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position_2px(position), phi_2px);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position_2py(position), phi_2py);
	EXPECT_FLOAT_EQ(hydrogenic.Wavefunction_Position_2pz(position), phi_2pz);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2s)
{
	// ARRANGE
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);

	std::function<double(libphysica::Vector)> integrand = [&hydrogenic](libphysica::Vector rVec) {
		Eigen::Vector3d r_aux(rVec[0], rVec[1], rVec[2]);
		double phi = hydrogenic.Wavefunction_Position_2s(r_aux);
		return phi * phi;
	};

	// ACT & ASSERT
	ASSERT_FLOAT_EQ(libphysica::Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2px)
{
	// ARRANGE
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	std::function<double(libphysica::Vector)> integrand = [&hydrogenic](libphysica::Vector rVec) {
		Eigen::Vector3d r_aux(rVec[0], rVec[1], rVec[2]);
		double phi = hydrogenic.Wavefunction_Position_2px(r_aux);
		return phi * phi;
	};

	// ACT & ASSERT
	ASSERT_FLOAT_EQ(libphysica::Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2py)
{
	// ARRANGE
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	std::function<double(libphysica::Vector)> integrand = [&hydrogenic](libphysica::Vector rVec) {
		Eigen::Vector3d r_aux(rVec[0], rVec[1], rVec[2]);
		double phi = hydrogenic.Wavefunction_Position_2py(r_aux);
		return phi * phi;
	};

	// ACT & ASSERT
	ASSERT_FLOAT_EQ(libphysica::Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI), 1.0);
}

TEST(TestHydrogenicWavefunctions, TestPositionWavefunctionNormalization2pz)
{
	// ARRANGE
	double Zeff = 2.0;
	Hydrogenic hydrogenic(Zeff, Zeff, Zeff);
	std::function<double(libphysica::Vector)> integrand = [&hydrogenic](libphysica::Vector rVec) {
		Eigen::Vector3d r_aux(rVec[0], rVec[1], rVec[2]);
		double phi = hydrogenic.Wavefunction_Position_2pz(r_aux);
		return phi * phi;
	};

	// ACT & ASSERT
	ASSERT_FLOAT_EQ(libphysica::Integrate_3D(integrand, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI), 1.0);
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
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum_2s(momentum)), phi_2s);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum_2px(momentum)), phi_2px);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum_2py(momentum)), phi_2py);
	EXPECT_FLOAT_EQ(std::abs(hydrogenic.Wavefunction_Momentum_2pz(momentum)), phi_2pz);
}