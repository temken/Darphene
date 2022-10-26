#include "graphene/Carbon_Wavefunctions.hpp"

#include <cmath>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"

namespace graphene
{
using namespace libphysica::natural_units;

// 1. Base class for atomic wavefunctions of carbon
double Carbon_Wavefunctions::Normalization_Position(const std::string& orbital)
{
	std::function<double(libphysica::Vector)> squared_wavefunction = [this, &orbital](libphysica::Vector rVec) {
		Eigen::Vector3d r_aux(rVec[0], rVec[1], rVec[2]);
		double phi = Wavefunction_Position(r_aux, orbital);
		return phi * phi;
	};
	return libphysica::Integrate_3D(squared_wavefunction, 0.0, 20.0 * Bohr_Radius, -1.0, 1.0, 0.0, 2.0 * M_PI);
}

double Carbon_Wavefunctions::Normalization_Momentum(const std::string& orbital)
{
	std::function<double(libphysica::Vector)> squared_wavefunction = [this, &orbital](libphysica::Vector kVec) {
		Eigen::Vector3d k_aux(kVec[0], kVec[1], kVec[2]);
		std::complex<double> phi = Wavefunction_Momentum(k_aux, orbital);
		return std::norm(phi);
	};
	return std::pow(2.0 * M_PI, -3.0) * libphysica::Integrate_3D(squared_wavefunction, 0.0, 1000.0 * keV, -1.0, 1.0, 0.0, 2.0 * M_PI, "Gauss-Kronrod");
}

// 2. Hydrogenic wavefunctions
// Position space wavefunctions
double Hydrogenic::Wavefunction_Position(const Eigen::Vector3d& xVec, const std::string& orbital)
{
	double r = xVec.norm();
	if(orbital == "2s")
	{
		double normalization = sqrt(Zeff_2s * Zeff_2s * Zeff_2s / 56.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * (1.0 - Zeff_2s * r / Bohr_Radius) * exp(-Zeff_2s * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2px")
	{
		double cos_theta	 = xVec[2] / r;
		double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
		double phi			 = atan2(xVec[1], xVec[0]);
		double normalization = sqrt(Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * cos(phi) * exp(-Zeff_2pxpy * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2py")
	{
		double cos_theta	 = xVec[2] / r;
		double sin_theta	 = sqrt(1.0 - cos_theta * cos_theta);
		double phi			 = atan2(xVec[1], xVec[0]);
		double normalization = sqrt(Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * sin_theta * sin(phi) * exp(-Zeff_2pxpy * r / 2.0 / Bohr_Radius);
	}
	else if(orbital == "2pz")
	{
		double cos_theta	 = xVec[2] / r;
		double normalization = sqrt(Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz / 32.0 / M_PI);
		return normalization / pow(Bohr_Radius, 1.5) * r / Bohr_Radius * cos_theta * exp(-Zeff_2pz * r / 2.0 / Bohr_Radius);
	}
	else
	{
		std::cerr << "Error in Hydrogenic::Wavefunction_Position(): Unknown orbital " << orbital << std::endl;
		exit(EXIT_FAILURE);
	}
}

// Momentum space wavefunctions
std::complex<double> Hydrogenic::Wavefunction_Momentum(const Eigen::Vector3d& kVec, const std::string& orbital)
{
	double k = kVec.norm();
	if(orbital == "2s")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff_2s * Zeff_2s * Zeff_2s * Zeff_2s * Zeff_2s);
		return normalization * pow(Bohr_Radius, 1.5) * (Bohr_Radius * Bohr_Radius * k * k - Zeff_2s * Zeff_2s / 4.0) / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2s * Zeff_2s / 4.0, 3.0);
	}
	else if(orbital == "2px")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[0] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pxpy * Zeff_2pxpy / 4.0, 3.0);
	}
	else if(orbital == "2py")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy * Zeff_2pxpy);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[1] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pxpy * Zeff_2pxpy / 4.0, 3.0);
	}
	else if(orbital == "2pz")
	{
		double normalization = sqrt(8.0 * M_PI * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz * Zeff_2pz);
		return normalization * pow(Bohr_Radius, 1.5) * Bohr_Radius * kVec[2] / pow(Bohr_Radius * Bohr_Radius * k * k + Zeff_2pz * Zeff_2pz / 4.0, 3.0);
	}
	else
	{
		std::cerr << "Error in Hydrogenic::Wavefunction_Position(): Unknown orbital " << orbital << std::endl;
		exit(EXIT_FAILURE);
	}
}

}	// namespace graphene
