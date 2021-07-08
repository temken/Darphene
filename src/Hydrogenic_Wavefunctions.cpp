#include "Hydrogenic_Wavefunctions.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

namespace graphene
{
using namespace libphysica::natural_units;

// Position space
double Hydrogenic_Wavefunction(const Eigen::Vector3d& position, const std::string& orbital, double Zeff)
{
	double r = position.norm();
	if(orbital == "2s")
	{
		double normalization = sqrt(Zeff * Zeff * Zeff / 56.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * (1.0 - Zeff * r / Bohr_Radius) * exp(-Zeff * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2px")
	{
		double cos_theta	 = position[2] / r;
		double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
		double phi			 = atan2(position[1], position[0]);
		double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * cos(phi) * exp(-Zeff * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2py")
	{
		double cos_theta	 = position[2] / r;
		double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
		double phi			 = atan2(position[1], position[0]);
		double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * sin(phi) * exp(-Zeff * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2pz")
	{
		double cos_theta	 = position[2] / r;
		double normalization = sqrt(Zeff * Zeff * Zeff * Zeff * Zeff / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * cos_theta * exp(-Zeff * r / 2.0 / Bohr_Radius);
	}
	else
	{
		std::cerr << "Error in Hydrogenic_Wavefunction(): Orbital " << orbital << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Momentum space
double Hydrogenic_Wavefunction_Momentum(const Eigen::Vector3d& momentum, const std::string& orbital, double Zeff)
{
	double k = momentum.norm();
	if(orbital == "2s")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff);
		return normalization * pow(Bohr_Radius, 1.5) * (Bohr_Radius * Bohr_Radius * k * k - Zeff * Zeff / 4.0) / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
	}
	else if(orbital == "2px")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[0] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
	}
	else if(orbital == "2py")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[1] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
	}
	else if(orbital == "2pz")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff * Zeff);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * momentum[2] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff * Zeff / 4.0, 3.0);
	}
	else
	{
		std::cerr << "Error in Hydrogenic_Wavefunction_Momentum(): Orbital " << orbital << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace graphene
