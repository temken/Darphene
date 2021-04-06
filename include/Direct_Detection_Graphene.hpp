#ifndef __Direct_Detection_Graphene_hpp_
#define __Direct_Detection_Graphene_hpp_

#include "Graphene.hpp"

// Headers form obscura
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"

namespace graphene
{
extern double perform_integral(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band);
extern double perform_integral_vegas(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band);

}	// namespace graphene

#endif