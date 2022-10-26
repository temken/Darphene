#include "graphene/Carbon_Wavefunctions.hpp"

#include <cmath>
#include <fstream>

#include "acb_hypgeom.h"

#include <boost/math/special_functions/factorials.hpp>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"

#include "version.hpp"

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

// 3. Roothaan-Hartree-Fock wavefunctions
double a0 = Bohr_Radius;
void Roothaan_Hartree_Fock::Import_RHF_Coefficients()
{
	std::vector<std::string> orbitals = {"2s", "2p"};
	for(int i = 0; i < 2; i++)
	{
		std::string filepath = TOP_LEVEL_DIR "data/C_" + orbitals[i] + ".txt";
		std::ifstream f;
		f.open(filepath);
		if(f.is_open())
		{
			double binding_energy;
			f >> binding_energy;
			double C, Z;
			unsigned int nin;
			while(f >> nin >> Z >> C)
			{
				n_lj[i].push_back(nin);
				Z_lj[i].push_back(Z);
				C_nlj[i].push_back(C);
			}
			f.close();
		}
	}
}

double Roothaan_Hartree_Fock::Radial_Wavefunction_Position(double r, const std::string& orbital) const
{
	int i		= orbital == "2s" ? 0 : 1;
	double R_nl = 0.0;
	for(unsigned int j = 0; j < C_nlj[i].size(); j++)
		R_nl += C_nlj[i][j] * std::pow(2.0 * Z_lj[i][j], n_lj[i][j] + 0.5) / sqrt(boost::math::factorial<double>(2.0 * n_lj[i][j])) * std::pow(r / a0, n_lj[i][j] - 1.0) * std::exp(-Z_lj[i][j] * r / a0);

	return std::pow(a0, -1.5) * R_nl;
}

std::complex<double> Hypergeometric_2F1(double a, double b, double c, double z)
{
	using namespace std::complex_literals;
	std::complex<double> result;
	slong prec;
	acb_t F_arb, a_arb, b_arb, c_arb, z_arb;
	acb_init(F_arb);
	acb_init(a_arb);
	acb_init(b_arb);
	acb_init(c_arb);
	acb_init(z_arb);
	acb_set_d(a_arb, a);
	acb_set_d(b_arb, b);
	acb_set_d(c_arb, c);
	acb_set_d(z_arb, z);
	for(prec = 80;; prec *= 2)
	{
		acb_hypgeom_2f1(F_arb, a_arb, b_arb, c_arb, z_arb, false, prec);
		if(acb_rel_accuracy_bits(F_arb) >= 53)
		{
			arb_t F_real, F_imag;
			arb_init(F_real);
			arb_init(F_imag);
			acb_get_real(F_real, F_arb);
			acb_get_imag(F_imag, F_arb);
			result = arf_get_d(arb_midref(F_real), ARF_RND_NEAR) + 1i * arf_get_d(arb_midref(F_imag), ARF_RND_NEAR);
			break;
		}
		else if(prec > 10000)
		{
			result = NAN;
			break;
		}
	}
	acb_clear(F_arb);
	acb_clear(a_arb);
	acb_clear(b_arb);
	acb_clear(c_arb);
	acb_clear(z_arb);
	return result;
}

std::complex<double> Roothaan_Hartree_Fock::Radial_Wavefunction_Momentum(double p, const std::string& orbital) const
{
	using namespace std::complex_literals;
	int i					  = orbital == "2s" ? 0 : 1;
	int l					  = orbital == "2s" ? 0 : 1;
	std::complex<double> R_nl = 0.0;
	for(unsigned int j = 0; j < C_nlj[i].size(); j++)
	{
		// Hypergeometric function 2F1(a,b,c,z)
		double a								= 0.5 * (2.0 + l + n_lj[i][j]);
		double b								= 0.5 * (3.0 + l + n_lj[i][j]);
		double c								= 1.5 + l;
		double z								= -std::pow(p * Bohr_Radius / Z_lj[i][j], 2.0);
		std::complex<double> hypergeometric_2F1 = Hypergeometric_2F1(a, b, c, z);

		R_nl += std::pow(2.0 * M_PI * Bohr_Radius / Z_lj[i][j], 1.5) * C_nlj[i][j] * std::pow(2.0, -l + n_lj[i][j]) * boost::math::factorial<double>(1.0 + n_lj[i][j] + l) / sqrt(boost::math::factorial<double>(2.0 * n_lj[i][j])) * std::pow(1i * p * Bohr_Radius / Z_lj[i][j], l) * hypergeometric_2F1 / tgamma(1.5 + l);
	}
	return R_nl;
}

Roothaan_Hartree_Fock::Roothaan_Hartree_Fock()
: C_nlj({{}, {}}), Z_lj({{}, {}}), n_lj({{}, {}})
{
	Import_RHF_Coefficients();
}

double Roothaan_Hartree_Fock::Wavefunction_Position(const Eigen::Vector3d& xVec, const std::string& orbital)
{
	double r = xVec.norm();
	if(orbital == "2s")
		return 1.0 / sqrt(4.0 * M_PI) * Radial_Wavefunction_Position(r, orbital);
	else if(orbital == "2px")
		return sqrt(3.0 / 4.0 / M_PI) * xVec[0] / r * Radial_Wavefunction_Position(r, orbital);
	else if(orbital == "2py")
		return sqrt(3.0 / 4.0 / M_PI) * xVec[1] / r * Radial_Wavefunction_Position(r, orbital);
	else if(orbital == "2pz")
		return sqrt(3.0 / 4.0 / M_PI) * xVec[2] / r * Radial_Wavefunction_Position(r, orbital);
	else
	{
		std::cout << "Error in Roothaan_Hartree_Fock::Wavefunction_Position(): Orbital " << orbital << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

std::complex<double> Roothaan_Hartree_Fock::Wavefunction_Momentum(const Eigen::Vector3d& kVec, const std::string& orbital)
{
	double k = kVec.norm();
	if(orbital == "2s")
		return 1.0 / sqrt(4.0 * M_PI) * Radial_Wavefunction_Momentum(kVec.norm(), orbital);
	else if(orbital == "2px")
		return sqrt(3.0 / 4.0 / M_PI) * kVec[0] / k * Radial_Wavefunction_Momentum(k, orbital);
	else if(orbital == "2py")
		return sqrt(3.0 / 4.0 / M_PI) * kVec[1] / k * Radial_Wavefunction_Momentum(k, orbital);
	else if(orbital == "2pz")
		return sqrt(3.0 / 4.0 / M_PI) * kVec[2] / k * Radial_Wavefunction_Momentum(k, orbital);
	else
	{
		std::cout << "Error in Roothaan_Hartree_Fock::Wavefunction_Momentum(): Orbital " << orbital << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace graphene
