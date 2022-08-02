#ifndef __Direct_Detection_Graphene_hpp_
#define __Direct_Detection_Graphene_hpp_

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle.hpp"

#include "graphene/DM_Particle_NREFT.hpp"
#include "graphene/Graphene.hpp"

namespace graphene
{
extern Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi);

// 1. Total scattering rates
// 1.1 Total rate per band
extern double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double R_Total_Standard_Simple(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
// 1.2 Total rate for all bands
extern double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points);
extern double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

// 2. Differential rate dR/dlnE
// 2.1 For one band
extern double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dlnE_Standard_Simple(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dlnE_NREFT(double Ee, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
// 2.2 For all bands
extern double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points);
extern double dR_dlnE_NREFT(double Ee, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

// 3. Differential rate dR / dcos(theta)
// 3.1 For one band
extern double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dcos_NREFT(double cos_theta, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);

// 3.2 For all bands
extern double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);
extern double dR_dcos_NREFT(double cos_theta, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

// 4. Differential rate d^2R/(dcos dphi)
// 4.1 For one band
extern double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dcos_dphi_NREFT(double cos_theta, double phi, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);

// 4.2 For all bands
extern double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);
extern double dR_dcos_dphi_NREFT(double cos_theta, double phi, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

// 5. Tabulation functions
// 5.1 Energy spectrum
extern std::vector<std::vector<double>> Tabulate_dR_dlnE_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points, int threads = 1);
extern std::vector<std::vector<double>> Tabulate_dR_dlnE_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads = 1);

// 5.2 Directional spectrum
extern std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads = 1);
extern std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads = 1);

// 5.3 Daily Modulation
extern libphysica::Vector Earth_Velocity(double t, double v_earth);
extern std::vector<std::vector<double>> Daily_Modulation_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads = 1);
extern std::vector<std::vector<double>> Daily_Modulation_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads = 1);

}	// namespace graphene

#endif