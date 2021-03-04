#ifndef __Hydrogenic_Wavefunctions_hpp_
#define __Hydrogenic_Wavefunctions_hpp_

#include <Eigen/Geometry>

namespace graphene
{

// Position space
extern double Hydrogenic_Wavefunction_2s(const Eigen::Vector3d& position, double Zeff);
extern double Hydrogenic_Wavefunction_2px(const Eigen::Vector3d& position, double Zeff);
extern double Hydrogenic_Wavefunction_2py(const Eigen::Vector3d& position, double Zeff);
extern double Hydrogenic_Wavefunction_2pz(const Eigen::Vector3d& position, double Zeff);

// Momentum space
extern double Hydrogenic_Wavefunction_Momentum_2s(const Eigen::Vector3d& momentum, double Zeff);
extern double Hydrogenic_Wavefunction_Momentum_2px(const Eigen::Vector3d& momentum, double Zeff);
extern double Hydrogenic_Wavefunction_Momentum_2py(const Eigen::Vector3d& momentum, double Zeff);
extern double Hydrogenic_Wavefunction_Momentum_2pz(const Eigen::Vector3d& momentum, double Zeff);

}	// namespace graphene

#endif