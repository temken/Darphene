#include "gtest/gtest.h"

#include <functional>
#include <limits>

#include "libphysica/Integration.hpp"
#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Natural_Units.hpp"

#include "Darphene/Carbon_Wavefunctions.hpp"

using namespace libphysica::natural_units;
using namespace Darphene;

// 2. Hydrogenic wavefunctions
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

// 3. Roothaan-Hartree-Fock wavefunctions
TEST(TestRoothaanHartreeFock, TestPositionWavefunctionValues)
{
	Eigen::Vector3d position({Bohr_Radius, Bohr_Radius / 2, Bohr_Radius / 3});
	Roothaan_Hartree_Fock rhf;
	double phi_2s  = -1.414450631e-9;
	double phi_2px = 1.877396747e-09;
	double phi_2py = 9.386983734e-10;
	double phi_2pz = 6.257989156e-10;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(rhf.Wavefunction_Position(position, "2s"), phi_2s);
	EXPECT_FLOAT_EQ(rhf.Wavefunction_Position(position, "2px"), phi_2px);
	EXPECT_FLOAT_EQ(rhf.Wavefunction_Position(position, "2py"), phi_2py);
	EXPECT_FLOAT_EQ(rhf.Wavefunction_Position(position, "2pz"), phi_2pz);
}

TEST(TestRoothaanHartreeFock, TestPositionWavefunctionNormalization)
{
	// ARRANGE
	double tol = 1e-4;
	Roothaan_Hartree_Fock rhf;
	std::vector<std::string> orbitals = {"2s", "2px", "2py", "2pz"};
	// ACT & ASSERT
	for(auto orbital : orbitals)
		EXPECT_NEAR(rhf.Normalization_Position(orbital), 1.0, tol);
}

TEST(TestRoothaanHartreeFock, TestMomentumWavefunctionValues)
{
	// ARRANGE
	Eigen::Vector3d momentum({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3 / Bohr_Radius});
	Roothaan_Hartree_Fock rhf;
	double phi_2s  = 2.354591573e7;
	double phi_2px = 1.0460452717e7;
	double phi_2py = 2.092090543e7;
	double phi_2pz = 3.138135815e7;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(std::abs(rhf.Wavefunction_Momentum(momentum, "2s")), phi_2s);
	EXPECT_FLOAT_EQ(std::abs(rhf.Wavefunction_Momentum(momentum, "2px")), phi_2px);
	EXPECT_FLOAT_EQ(std::abs(rhf.Wavefunction_Momentum(momentum, "2py")), phi_2py);
	EXPECT_FLOAT_EQ(std::abs(rhf.Wavefunction_Momentum(momentum, "2pz")), phi_2pz);
}

TEST(TestRoothaanHartreeFock, TestMomentumWavefunctionNormalization)
{
	// ARRANGE
	double tol = 1e-4;
	Roothaan_Hartree_Fock rhf;
	std::vector<std::string> orbitals = {"2s", "2px", "2py", "2pz"};
	// ACT & ASSERT
	for(auto orbital : orbitals)
		EXPECT_NEAR(rhf.Normalization_Momentum(orbital), 1.0, tol);
}