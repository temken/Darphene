#ifndef __Direct_Detection_Graphene_hpp_
#define __Direct_Detection_Graphene_hpp_

#include "Graphene.hpp"

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"

namespace graphene
{
extern Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi);
extern libphysica::Vector Earth_Velocity(double t, double v_earth);

// Correct response functions, simplified velocity integral
extern double dR_dlnE_simplified(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total_simplified(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total_simplified(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method = "Vegas");

// Correct response functions, full velocity integral
extern double R_Total(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method = "Vegas");
extern double dR_dlnE(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double dR_dcosk_dphik(double cos_k, double phi_k, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");

extern std::vector<std::vector<double>> Tabulate_dR_dlnE(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method = "Full");
extern std::vector<std::vector<double>> Tabulate_dR_dcosk_dphik(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene);

}	// namespace graphene

#endif