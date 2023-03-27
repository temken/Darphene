#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "Darphene/DM_Particle_NREFT.hpp"

using namespace Darphene;
using namespace libphysica::natural_units;

TEST(TestDMFormFactor, TestDefaultConstructor)
{
	// ARRANGE
	double q = keV;

	// ACT
	DM_Form_Factor ff;

	// ASSERT
	ASSERT_FLOAT_EQ(ff(q), 1.0);
}

TEST(TestDMFormFactor, TestFormFactors)
{
	// ARRANGE
	double q								   = keV;
	double qRef								   = aEM * mElectron;
	double param							   = 3.0;
	std::vector<std::string> form_factor_types = {"Contact", "Long-Range", "General", "Power"};
	std::vector<double> results				   = {1.0, qRef * qRef / q / q, (qRef * qRef + param * param) / (q * q + param * param), std::pow(q / qRef, param)};

	// ACT && ASSERT
	std::vector<DM_Form_Factor> DM_Form_Factors = {};
	for(unsigned int i = 0; i < form_factor_types.size(); i++)
		DM_Form_Factors.push_back(DM_Form_Factor(form_factor_types[i], param));

	// ASSERT
	for(unsigned int i = 0; i < form_factor_types.size(); i++)
		EXPECT_FLOAT_EQ(DM_Form_Factors[i](q), results[i]);
}

TEST(TestDMFormFactor, TestPrintSummary)
{
	// ARRANGE
	DM_Form_Factor ff;

	// ACT & ASSERT
	ff.Print_Summary();
}

// 2. Class for the DM particle interacting with electrons via effective operators O_i
TEST(TestDMParticleNREFT, TestConstructors)
{
	// ARRANGE
	double mDM = 100.0 * MeV;
	double s   = 0.5;

	// ACT
	DM_Particle_NREFT DM(mDM, s);

	// ASSERT
	EXPECT_FLOAT_EQ(DM.mass, mDM);
	EXPECT_FLOAT_EQ(DM.spin, s);
}

TEST(TestDMParticleNREFT, SetCoupling)
{
	// ARRANGE
	double mDM = 100.0 * MeV;
	double s   = 0.5;
	DM_Particle_NREFT DM(mDM, s);
	double mu	 = libphysica::Reduced_Mass(mDM, mElectron);
	double c	 = 1e-3;
	double sigma = c * c * mu * mu / 16.0 / M_PI / mDM / mDM / mElectron / mElectron;
	// ACT
	DM.Set_Coupling(1, c, "Contact", 0.0);

	// ASSERT
	EXPECT_FLOAT_EQ(DM.Sigma_Electron(), sigma);
}

TEST(TestDMParticleNREFT, SetCrossSection)
{
	// ARRANGE
	double mDM = 100.0 * MeV;
	double s   = 0.5;
	DM_Particle_NREFT DM(mDM, s);
	double sigma = 1e-37 * cm * cm;
	// ACT
	DM.Set_Cross_Section(1, sigma, "Contact", 0.0);

	// ASSERT
	EXPECT_FLOAT_EQ(DM.Sigma_Electron(), sigma);
}

TEST(TestDMParticleNREFT, TestResponseFunction)
{
	// ARRANGE
	double mDM = 100.0 * MeV;
	double s   = 0.5;
	DM_Particle_NREFT DM(mDM, s);
	double c1 = 1e-3;
	DM.Set_Coupling(1, c1, "Contact", 0.0);

	Eigen::Vector3d qVec   = {1.0 * keV, 1.0 * keV, 1.0 * keV};
	Eigen::Vector3d velDM  = {0.0, 0.0, 1e-3};
	Eigen::Vector3d kPrime = {0.0, 0.0, 0.0};

	// ACT
	double R = DM.Response_Function(qVec, velDM, kPrime);

	// ASSERT
	EXPECT_FLOAT_EQ(R, c1 * c1);
}

TEST(TestDMParticleNREFT, PrintSummary)
{
	// ARRANGE
	double mDM = 100.0 * MeV;
	double s   = 0.5;
	DM_Particle_NREFT DM(mDM, s);

	// ACT & ASSERT
	DM.Print_Summary();
}

// 3. Class for benchmark models
TEST(TestDMParticleNREFT, TestDarkPhoton)
{
	// ARRANGE
	double mDM			 = 100.0 * MeV;
	double sigma_e		 = 1.0 * pb;
	DM_Particle_NREFT DM = DM_Dark_Photon(mDM, sigma_e);

	// ACT & ASSERT
	EXPECT_FLOAT_EQ(DM.mass, mDM);
	EXPECT_FLOAT_EQ(DM.Sigma_Electron(), sigma_e);
}