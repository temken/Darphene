#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "graphene/DM_Particle_NREFT.hpp"

using namespace graphene;
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
	for(int i = 0; i < form_factor_types.size(); i++)
		DM_Form_Factors.push_back(DM_Form_Factor(form_factor_types[i], param));

	// ASSERT
	for(int i = 0; i < form_factor_types.size(); i++)
		EXPECT_FLOAT_EQ(DM_Form_Factors[i](q), results[i]);
}

TEST(TestDMFormFactor, TestPrintSummary)
{
	// ARRANGE
	DM_Form_Factor ff;

	// ACT & ASSERT
	ff.Print_Summary();
}
