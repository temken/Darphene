#ifndef __Direct_Detection_Standard_hpp_
#define __Direct_Detection_Standard_hpp_

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle.hpp"

#include "graphene/Graphene.hpp"

namespace graphene
{
// 0. Auxiliary functions
extern Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi);
extern libphysica::Vector Earth_Velocity(double t, double v_earth);
extern double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum, double work_function);

// 1. Rates and spectra for a single band
// 1.1 Total rate per band
extern double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double R_Total_Standard_Simple(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);

// 1.2 Differential rates per band
extern double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dlnE_Standard_Simple(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);
extern double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points);

// 2. Rates and spectra for all bands
// 2.1 Total rate
extern double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points);

// 2.2 Differential rates
extern double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points);
extern double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);
extern double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

// 3. Tabulation functions
extern std::vector<std::vector<double>> Tabulate_dR_dlnE_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points);
extern std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);
extern std::vector<std::vector<double>> Daily_Modulation_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

}	// namespace graphene

#endif