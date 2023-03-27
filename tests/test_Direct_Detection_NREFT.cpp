#include "gtest/gtest.h"

#include "Darphene/Direct_Detection_NREFT.hpp"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle_Standard.hpp"

#include "Darphene/Direct_Detection_Standard.hpp"

using namespace Darphene;
using namespace libphysica::natural_units;

// TEST(TestDirectDetectionNREFT, TestRTotalNREFT)
// {
// 	// ARRANGE
// 	// extern double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
// 	double mDM = 250.0 * MeV;
// 	DM_Particle_NREFT DM(mDM);
// 	DM.Set_Coupling(1, 1.0);
// 	obscura::Standard_Halo_Model SHM;
// 	Graphene graphene;
// 	int MC_points  = 1000;
// 	double rel_tol = 0.3;

// 	// ACT
// 	double R_pi	 = R_Total_NREFT(DM, SHM, graphene, 0, MC_points);
// 	double R_s1	 = R_Total_NREFT(DM, SHM, graphene, 1, MC_points);
// 	double R_s2	 = R_Total_NREFT(DM, SHM, graphene, 2, MC_points);
// 	double R_s3	 = R_Total_NREFT(DM, SHM, graphene, 3, MC_points);
// 	double R_tot = R_Total_NREFT(DM, SHM, graphene, MC_points);

// 	// ASSERT
// 	EXPECT_NEAR(R_pi + R_s1 + R_s2 + R_s3, R_tot, rel_tol * R_tot);
// }

// TEST(TestDirectDetectionNREFT, TestRTotalNREFTvsStandard)
// {
// 	// ARRANGE
// 	double mDM			 = 250.0 * MeV;
// 	double sigma		 = pb;
// 	DM_Particle_NREFT DM = DM_Dark_Photon(mDM, sigma);

// 	obscura::DM_Particle_SI DM_SI(mDM);
// 	DM_SI.Set_Sigma_Electron(sigma);

// 	obscura::Standard_Halo_Model SHM;
// 	Graphene graphene;
// 	int MC_points  = 1000;
// 	double rel_tol = 0.3;

// 	// ACT
// 	double R_pi_Standard = R_Total_Standard(DM_SI, SHM, graphene, 0, MC_points);
// 	double R_pi_NREFT	 = R_Total_NREFT(DM, SHM, graphene, 0, MC_points);

// 	// ASSERT
// 	EXPECT_NEAR(R_pi_Standard, R_pi_NREFT, rel_tol * R_pi_NREFT);
// }