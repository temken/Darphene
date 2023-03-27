#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "graphene/Configuration.hpp"

using namespace Darphene;
using namespace libphysica::natural_units;

TEST(TestConfiguration, TestConstructorNREFT)
{
	// ARRANGE
	// ACT
	Configuration cfg(".test_config_NREFT.cfg");

	// ASSERT
	EXPECT_EQ(cfg.run_modus, "Custom");
	EXPECT_EQ(cfg.carbon_wavefunctions, "RHF");
	EXPECT_EQ(cfg.MC_points, 1000);
	EXPECT_EQ(cfg.grid_points, 100);
	EXPECT_FLOAT_EQ(cfg.graphene_work_function, 4.3 * eV);
	EXPECT_FLOAT_EQ(cfg.time, 6.0 * hr);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->mass, 0.1 * GeV);
}

TEST(TestConfiguration, TestConstructorAnapole)
{
	// ARRANGE
	// ACT
	Configuration cfg(".test_config_Anapole.cfg");

	// ASSERT
	EXPECT_EQ(cfg.run_modus, "Custom");
	EXPECT_EQ(cfg.carbon_wavefunctions, "RHF");
	EXPECT_EQ(cfg.MC_points, 1000);
	EXPECT_EQ(cfg.grid_points, 100);
	EXPECT_FLOAT_EQ(cfg.graphene_work_function, 4.3 * eV);
	EXPECT_FLOAT_EQ(cfg.time, 6.0 * hr);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->mass, 0.1 * GeV);
}

TEST(TestConfiguration, TestConstructorDarkPhoton)
{
	// ARRANGE
	// ACT
	Configuration cfg(".test_config_Dark_Photon.cfg");

	// ASSERT
	EXPECT_EQ(cfg.run_modus, "Custom");
	EXPECT_EQ(cfg.carbon_wavefunctions, "RHF");
	EXPECT_EQ(cfg.MC_points, 1000);
	EXPECT_EQ(cfg.grid_points, 100);
	EXPECT_FLOAT_EQ(cfg.graphene_work_function, 4.3 * eV);
	EXPECT_FLOAT_EQ(cfg.time, 6.0 * hr);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->mass, 0.1 * GeV);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->Sigma_Electron(), 1.0e-37 * cm * cm);
}

TEST(TestConfiguration, TestConstructorElectricDipole)
{
	// ARRANGE
	// ACT
	Configuration cfg(".test_config_Electric_Dipole.cfg");

	// ASSERT
	EXPECT_EQ(cfg.run_modus, "Custom");
	EXPECT_EQ(cfg.carbon_wavefunctions, "RHF");
	EXPECT_EQ(cfg.MC_points, 1000);
	EXPECT_EQ(cfg.grid_points, 100);
	EXPECT_FLOAT_EQ(cfg.graphene_work_function, 4.3 * eV);
	EXPECT_FLOAT_EQ(cfg.time, 6.0 * hr);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->mass, 0.1 * GeV);
}

TEST(TestConfiguration, TestConstructorMagneticDipole)
{
	// ARRANGE
	// ACT
	Configuration cfg(".test_config_Magnetic_Dipole.cfg");

	// ASSERT
	EXPECT_EQ(cfg.run_modus, "Custom");
	EXPECT_EQ(cfg.carbon_wavefunctions, "RHF");
	EXPECT_EQ(cfg.MC_points, 1000);
	EXPECT_EQ(cfg.grid_points, 100);
	EXPECT_FLOAT_EQ(cfg.graphene_work_function, 4.3 * eV);
	EXPECT_FLOAT_EQ(cfg.time, 6.0 * hr);
	EXPECT_FLOAT_EQ(cfg.DM_NREFT->mass, 0.1 * GeV);
}