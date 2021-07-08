#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

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

TEST(TestGraphene, TestEnergyDispersionSigma)
{
	//ARRANGE
	double aCC										  = 1.42 * Angstrom;
	double a										  = aCC * sqrt(3.0);
	std::vector<Eigen::Vector3d> high_symmetry_points = {{0.0, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0}, {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0}};
	Graphene graphene;
	std::vector<std::vector<double>> results = {
		{-17.833129584352082 * eV, -2.93125 * eV, -2.93125 * eV, 3.08466 * eV, 3.08466 * eV, 31.425824175824184 * eV},
		{-14.802126408956223 * eV, -11.691743504673768 * eV, -7.06817 * eV, 10.476314847203604 * eV, 12.661549197487779 * eV, 20.1867 * eV},
		{-13.064658876661905 * eV, -13.064658876661898 * eV, -8.56991 * eV, 12.119618979416404 * eV, 12.119618979416412 * eV, 20.604255319148937 * eV}};
	double tolerance = 1.0e-5 * eV;
	//ACT & ASSERT
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
	Graphene graphene;
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

TEST(TestGraphene, TestWaveFunctionPi)
{
	// ARRANGE
	Graphene graphene;

	Eigen::Vector3d rVec({Bohr_Radius, 2 * Bohr_Radius, 3 * Bohr_Radius});
	Eigen::Vector3d lVec({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3.0 / Bohr_Radius});

	double psiNorm = 6.1979866e-11;
	// ACT & ASSERT
	ASSERT_FLOAT_EQ(std::abs(graphene.Wavefunction_Pi(rVec, lVec, 0)), psiNorm);
}

TEST(TestGraphene, TestWaveFunctionPiAnalytic)
{
	// ARRANGE
	Graphene graphene;

	Eigen::Vector3d rVec({Bohr_Radius, 2 * Bohr_Radius, 3 * Bohr_Radius});
	Eigen::Vector3d lVec({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3.0 / Bohr_Radius});

	double psi			= std::abs(graphene.Wavefunction_Pi(rVec, lVec, 0));
	double psi_analytic = std::abs(graphene.Wavefunction_Pi_Analytic(rVec, lVec));
	// ACT & ASSERT
	ASSERT_FLOAT_EQ(psi, psi_analytic);
}

TEST(TestGraphene, TestWaveFunctionSigma)
{
	// ARRANGE
	Graphene graphene;

	Eigen::Vector3d rVec({Bohr_Radius, 2 * Bohr_Radius, 3 * Bohr_Radius});
	Eigen::Vector3d lVec({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3.0 / Bohr_Radius});

	std::vector<double> psiNorm = {2.41481e-11, 6.1055726e-12, 2.8370941e-12, 5.7211831e-12,
								   7.28456e-12, 1.0713087e-11};
	// ACT & ASSERT
	for(int i = 0; i < 6; i++)
		EXPECT_FLOAT_EQ(std::abs(graphene.Wavefunction_Sigma(rVec, lVec, i)), psiNorm[i]);
}