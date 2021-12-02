#ifndef __Direct_Detection_Graphene_hpp_
#define __Direct_Detection_Graphene_hpp_

#include "Graphene.hpp"

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"

namespace graphene
{
extern double dR_dlnE(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method = "Vegas");

}	// namespace graphene

#endif