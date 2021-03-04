#include "Hydrogenic_Wavefunctions.hpp"

#include <cmath>

//Headers from libphysica
#include "Natural_Units.hpp"

namespace graphene
{
using namespace libphysica::natural_units;

// Position space
double Hydrogenic_Wavefunction_2s(const Eigen::Vector3d& position, double Zeff)
{
	double r = position.norm();

	double normalization = sqrt(Zeff * Zeff * Zeff / 56.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * (1.0 - Zeff * r / Bohr_Radius) * exp(-Zeff * r / 2.0 / Bohr_Radius);
}

double Hydrogenic_Wavefunction_2px(const Eigen::Vector3d& position, double Zeff)
{
	double r		 = position.norm();
	double cos_theta = position[2] / r;
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	double phi		 = atan2(position[1], position[0]);

	double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * cos(phi) * exp(-Zeff * r / 2.0 / Bohr_Radius);
}

double Hydrogenic_Wavefunction_2py(const Eigen::Vector3d& position, double Zeff)
{
	double r		 = position.norm();
	double cos_theta = position[2] / r;
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	double phi		 = atan2(position[1], position[0]);

	double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * sin(phi) * exp(-Zeff * r / 2.0 / Bohr_Radius);
}

double Hydrogenic_Wavefunction_2pz(const Eigen::Vector3d& position, double Zeff)
{
	double r		 = position.norm();
	double cos_theta = position[2] / r;

	double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
	return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * cos_theta * exp(-Zeff * r / 2.0 / Bohr_Radius);
}

// Momentum space
double Hydrogenic_Wavefunction_Momentum_2s(const Eigen::Vector3d& momentum, double Zeff)
{
	double k = momentum.norm();

	double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff);
	return normalization * pow(Bohr_Radius, 1.5) * (Bohr_Radius * Bohr_Radius * k * k - Zeff * Zeff / 4.0) / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
}

double Hydrogenic_Wavefunction_Momentum_2px(const Eigen::Vector3d& momentum, double Zeff)
{
	double k = momentum.norm();

	double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[0] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
}

double Hydrogenic_Wavefunction_Momentum_2py(const Eigen::Vector3d& momentum, double Zeff)
{
	double k = momentum.norm();

	double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[1] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
}

double Hydrogenic_Wavefunction_Momentum_2pz(const Eigen::Vector3d& momentum, double Zeff)
{
	double k = momentum.norm();

	double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
	return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[2] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
}

}	// namespace graphene
