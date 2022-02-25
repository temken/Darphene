#ifndef __Direct_Detection_Graphene_hpp_
#define __Direct_Detection_Graphene_hpp_

#include "Graphene.hpp"

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"

namespace graphene
{
extern Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi);
extern libphysica::Vector Earth_Velocity(double t, double v_earth);

// Reproducing Hochberg et al.
extern double dR_dlnE(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");

// Correct response functions, simplified velocity integral
extern double dR_dlnE_corrected(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total_corrected(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total_corrected(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method = "Vegas");

// Correct response functions, full velocity integral
extern double R_Total_Full_Integral(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");
extern double R_Total_Full_Integral(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method = "Vegas");
extern double dR_dlnE_Full_Integral(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");

}	// namespace graphene

#endif