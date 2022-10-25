#include "graphene/Carbon_Wavefunctions.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

namespace graphene
{
using namespace libphysica::natural_units;

// 2. Hydrogenic wavefunctions
// Position space wavefunctions
double Hydrogenic::Wavefunction_Position_2s(const Eigen::Vector3d& xVec)
{
	double r			 = xVec.norm();
	double normalization = sqrt(Zeff_2s * Zeff_2s * Zeff_2s / 56.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * (1.0 - Zeff_2s * r / Bohr_Radius) * exp(-Zeff_2s * r / 2.0 / Bohr_Radius);
}
double Hydrogenic::Wavefunction_Position_2px(const Eigen::Vector3d& xVec)
{
	double r			 = xVec.norm();
	double cos_theta	 = xVec[2] / r;
	double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
	double phi			 = atan2(xVec[1], xVec[0]);
	double normalization = sqrt(Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * cos(phi) * exp(-Zeff_2pxpy * r / 2.0 / Bohr_Radius);
}
double Hydrogenic::Wavefunction_Position_2py(const Eigen::Vector3d& xVec)
{
	double r			 = xVec.norm();
	double cos_theta	 = xVec[2] / r;
	double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
	double phi			 = atan2(xVec[1], xVec[0]);
	double normalization = sqrt(Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * sin(phi) * exp(-Zeff_2pxpy * r / 2.0 / Bohr_Radius);
}
double Hydrogenic::Wavefunction_Position_2pz(const Eigen::Vector3d& xVec)
{
	double r			 = xVec.norm();
	double cos_theta	 = xVec[2] / r;
	double normalization = sqrt(Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * cos_theta * exp(-Zeff_2pz * r / 2.0 / Bohr_Radius);
}

// Momentum space wavefunctions
std::complex<double> Hydrogenic::Wavefunction_Momentum_2s(const Eigen::Vector3d& kVec)
{
	double k			 = kVec.norm();
	double normalization = sqrt(8.0 * M_PI * Zeff_2s * Zeff_2s * Zeff_2s * Zeff_2s * Zeff_2s);
	return normalization * pow(Bohr_Radius, 1.5) * (Bohr_Radius * Bohr_Radius * k * k - Zeff_2s * Zeff_2s / 4.0) / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2s * Zeff_2s / 4.0, 3.0);
}
std::complex<double> Hydrogenic::Wavefunction_Momentum_2px(const Eigen::Vector3d& kVec)
{
	double k			 = kVec.norm();
	double normalization = sqrt(8.0 * M_PI * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[0] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pxpy * Zeff_2pxpy / 4.0, 3.0);
}
std::complex<double> Hydrogenic::Wavefunction_Momentum_2py(const Eigen::Vector3d& kVec)
{
	double k			 = kVec.norm();
	double normalization = sqrt(8.0 * M_PI * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[1] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pxpy * Zeff_2pxpy / 4.0, 3.0);
}
std::complex<double> Hydrogenic::Wavefunction_Momentum_2pz(const Eigen::Vector3d& kVec)
{
	double k			 = kVec.norm();
	double normalization = sqrt(8.0 * M_PI * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[2] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pz * Zeff_2pz / 4.0, 3.0);
}

}	// namespace graphene
