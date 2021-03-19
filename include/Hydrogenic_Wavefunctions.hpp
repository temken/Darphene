#ifndef __Hydrogenic_Wavefunctions_hpp_
#define __Hydrogenic_Wavefunctions_hpp_

#include <Eigen/Geometry>
#include <string>

namespace graphene
{

// Position space
extern double Hydrogenic_Wavefunction(const Eigen::Vector3d& position, const std::string& orbital, double Zeff);

// Momentum space
extern double Hydrogenic_Wavefunction_Momentum(const Eigen::Vector3d& momentum, const std::string& orbital, double Zeff);

}	// namespace graphene

#endif