#include "graphene/Direct_Detection_Standard.hpp"

#include <mpi.h>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"

namespace Darphene
{

using namespace libphysica::natural_units;

// 0. Auxiliary functions
Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi)
{
	if(r == 0)
		return Eigen::Vector3d(0, 0, 0);
	else
		return Eigen::Vector3d(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
}

Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi, const Eigen::Vector3d& axis)
{
	Eigen::Vector3d ev = axis.normalized();
	if(ev[2] == 1.0 || axis.norm() == 0.0)
		return Spherical_Coordinates(r, theta, phi);
	else
	{
		double aux = sqrt(1.0 - pow(ev[2], 2.0));

		double cos_theta = cos(theta);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
		double cos_phi	 = cos(phi);
		double sin_phi	 = sin(phi);

		Eigen::Vector3d unit_vector({cos_theta * ev[0] + sin_theta / aux * (ev[0] * ev[2] * cos_phi - ev[1] * sin_phi),
									 cos_theta * ev[1] + sin_theta / aux * (ev[1] * ev[2] * cos_phi + ev[0] * sin_phi),
									 cos_theta * ev[2] - aux * cos_phi * sin_theta});
		return r * unit_vector;
	}
}

libphysica::Vector Earth_Velocity(double t, double v_earth)
{
	double alpha = 42.0 * deg;
	double beta	 = 2.0 * M_PI * t / day;

	double cosa = cos(alpha);
	double sina = sin(alpha);
	double cosb = cos(beta);
	double sinb = sin(beta);

	return libphysica::Vector({v_earth * sina * sinb, v_earth * sina * cosa * (cosb - 1.0), v_earth * (cosa * cosa + sina * sina * cosb)});
}

double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum, double work_function)
{
	double E_final = final_momentum * final_momentum / 2.0 / mElectron;
	return (E_final - energy_crystal + work_function) / q + q / 2.0 / mDM;
}

// 1. Rates and spectra for a single band
// 1.1 Total rate per band
double R_Total_Standard(const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = 1.0 / (2.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double kf			= x[3];
		double cos_theta_kf = x[4];
		double phi_kf		= x[5];
		double cos_theta_v	= x[6];
		double phi_v		= x[7];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);

		// Determine the angle between vVec and qVec
		Eigen::Vector3d v_unitvector = Spherical_Coordinates(1.0, acos(cos_theta_v), phi_v);
		double cos_alpha			 = v_unitvector.dot(qVec) / q;

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, cos_theta_kf_min, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, cos_theta_kf_max, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");

	return std::isnan(result) ? 0.0 : result;
}

double R_Total_Standard_Simple(const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = 1.0 / (4.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, band, mDM](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double kf			= x[3];
		double cos_theta_kf = x[4];
		double phi_kf		= x[5];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k	= graphene.Valence_Band_Energies(kVec, band);
		double vMin = vMinimum_Graphene(mDM, q, E_k, kf, graphene.work_function);

		return kf * kf * q * DM_distr.Eta_Function(vMin) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, cos_theta_kf_min, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, cos_theta_kf_max, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 1.2 Differential rates per band
double dR_dlnE_Standard(double Ee, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = std::sqrt(2.0 * std::pow(mElectron * Ee, 3.0)) / (2.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;
	double kf		 = std::sqrt(2.0 * mElectron * Ee);

	double vMax = DM_distr.Maximum_DM_Speed();
	double EMax = mDM / 2.0 * vMax * vMax;
	if(Ee > EMax)
		return 0.0;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, kf, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double cos_theta_kf = x[3];
		double phi_kf		= x[4];
		double cos_theta_v	= x[5];
		double phi_v		= x[6];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);

		// Determine the angle between vVec and qVec
		Eigen::Vector3d v_unitvector = Spherical_Coordinates(1.0, acos(cos_theta_v), phi_v);
		double cos_alpha			 = v_unitvector.dot(qVec) / q;

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, cos_theta_kf_min, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, cos_theta_kf_max, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

double dR_dlnE_Standard_Simple(double Ee, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = std::sqrt(2.0 * std::pow(mElectron * Ee, 3.0)) / (4.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;
	double kf		 = std::sqrt(2.0 * mElectron * Ee);

	double vMax = DM_distr.Maximum_DM_Speed();
	double EMax = mDM / 2.0 * vMax * vMax;
	if(Ee > EMax)
		return 0.0;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, band, mDM, kf](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double cos_theta_kf = x[3];
		double phi_kf		= x[4];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k	= graphene.Valence_Band_Energies(kVec, band);
		double vMin = vMinimum_Graphene(mDM, q, E_k, kf, graphene.work_function);

		return q * DM_distr.Eta_Function(vMin) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, cos_theta_kf_min, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, cos_theta_kf_max, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

double dR_dcos_Standard(double cos_theta, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = 1.0 / (2.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [cos_theta, &DM_distr, &graphene, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
		double q		   = x[0];
		double cos_theta_q = x[1];
		double phi_q	   = x[2];
		double kf		   = x[3];
		double phi_kf	   = x[4];
		double cos_theta_v = x[5];
		double phi_v	   = x[6];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta), phi_kf);

		// Determine the angle between vVec and qVec
		Eigen::Vector3d v_unitvector = Spherical_Coordinates(1.0, acos(cos_theta_v), phi_v);
		double cos_alpha			 = v_unitvector.dot(qVec) / q;

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");

	return std::isnan(result) ? 0.0 : result;
}

double dR_dcos_dphi_Standard(double cos_theta, double phi, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double prefactor = 1.0 / (2.0 * M_PI) * DM_distr.DM_density / mDM * graphene.N_cell * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, band, mDM, vMax, cos_theta, phi](const std::vector<double>& x, const double wgt) {
		double q		   = x[0];
		double cos_theta_q = x[1];
		double phi_q	   = x[2];
		double kf		   = x[3];
		double cos_theta_v = x[4];
		double phi_v	   = x[5];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta), phi);

		// Determine the angle between vVec and qVec
		Eigen::Vector3d v_unitvector = Spherical_Coordinates(1.0, acos(cos_theta_v), phi_v);
		double cos_alpha			 = v_unitvector.dot(qVec) / q;

		// Determine the crystal momentum vector k =  l^|| - G* (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
		Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
		Eigen::Vector3d kVec = l_parallel - G;

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 2. Rates and spectra for all bands
// 2.1 Total rate
double R_Total_Standard(const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points)
{
	double rTot = 0.0;
	for(int band = 0; band < 4; band++)
	{
		if(velocity_integral == "Full")
			rTot += R_Total_Standard(DM, DM_distr, graphene, band, MC_points);
		else if(velocity_integral == "Simplified")
			rTot += R_Total_Standard_Simple(DM, DM_distr, graphene, band, MC_points);
	}
	return rTot;
}

// 2.2 Differential rates
double dR_dlnE_Standard(double Ee, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points)
{
	double dRdlnE = 0.0;
	for(int band = 0; band < 4; band++)
	{
		if(velocity_integral == "Full")
			dRdlnE += dR_dlnE_Standard(Ee, DM, DM_distr, graphene, band, MC_points);
		else if(velocity_integral == "Simplified")
			dRdlnE += dR_dlnE_Standard_Simple(Ee, DM, DM_distr, graphene, band, MC_points);
	}
	return dRdlnE;
}

double dR_dcos_Standard(double cos_theta, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_Standard(cos_theta, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

double dR_dcos_dphi_Standard(double cos_theta, double phi, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_dphi_Standard(cos_theta, phi, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

// 3. Tabulation functions
std::vector<std::vector<double>> Tabulate_dR_dlnE_Standard(int points, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points)
{
	double E_min			   = 0.1 * eV;
	double v_max			   = DM_distr.Maximum_DM_Speed();
	double E_max			   = 0.99 * DM.mass / 2.0 * v_max * v_max;
	std::vector<double> E_kist = libphysica::Log_Space(E_min, E_max, points);
	std::vector<std::vector<double>> spectrum;
	for(auto& E_er : E_kist)
	{
		std::vector<double> row = {E_er};
		double dRdlnE			= 0.0;
		for(int band = 0; band < 4; band++)
		{
			double band_contribution = (velocity_integral == "Full") ? dR_dlnE_Standard(E_er, DM, DM_distr, graphene, band, MC_points) : dR_dlnE_Standard_Simple(E_er, DM, DM_distr, graphene, band, MC_points);
			row.push_back(band_contribution);
			dRdlnE += band_contribution;
		}
		row.push_back(dRdlnE);
		spectrum.push_back(row);
	}
	return spectrum;
}

std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_Standard(int points, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	auto cos_k_list = libphysica::Linear_Space(-1.0, 1.0, points);
	auto phi_k_list = libphysica::Linear_Space(0.0, 2.0 * M_PI, points);

	std::vector<std::vector<double>> spectrum;
	for(auto& cos_theta : cos_k_list)
		for(auto& phi : phi_k_list)
			spectrum.push_back({cos_theta, phi, dR_dcos_dphi_Standard(cos_theta, phi, DM, DM_distr, graphene, MC_points)});

	return spectrum;
}

std::vector<std::vector<double>> Daily_Modulation_Standard(int points, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	// Total rate over the course of a day
	std::vector<double> t_list = libphysica::Linear_Space(0.0, 24.0, points);
	std::vector<std::vector<double>> daily_modulation_list;
	double vEarth = dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Get_Observer_Velocity().Norm();
	for(auto& t : t_list)
	{
		dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Set_Observer_Velocity(Earth_Velocity(t * hr, vEarth));
		double R = R_Total_Standard(DM, DM_distr, graphene, "Full", MC_points);
		daily_modulation_list.push_back({t, R});
	}
	return daily_modulation_list;
}

// 4. DM detector class

DM_Detector_Graphene::DM_Detector_Graphene(std::string& label, double exposure)
: obscura::DM_Detector(label, exposure, "Electrons")
{
}

double DM_Detector_Graphene::Maximum_Energy_Deposit(obscura::DM_Particle& DM, const obscura::DM_Distribution& DM_distr) const
{
	double vMax = DM_distr.Maximum_DM_Speed();
	return DM.mass * vMax * vMax / 2.0;
}

double DM_Detector_Graphene::Minimum_DM_Speed(obscura::DM_Particle& DM) const
{
	double threshold = graphene.work_function;
	double vMin		 = sqrt(2.0 * threshold / DM.mass);
	return vMin;
}

double DM_Detector_Graphene::Minimum_DM_Mass(obscura::DM_Particle& DM, const obscura::DM_Distribution& DM_distr) const
{
	double threshold = graphene.work_function;
	double vMax		 = DM_distr.Maximum_DM_Speed();
	double mMin		 = 2.0 * threshold / vMax / vMax;
	return mMin;
}

double DM_Detector_Graphene::dRdE(double E, const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr)
{
	int MC_points = 10000;
	return 1.0 / E * dR_dlnE_Standard(E, DM, DM_distr, graphene, "Full", MC_points);
}

double DM_Detector_Graphene::DM_Signals_Total(const obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr)
{
	int MC_points = 10000;
	return R_Total_Standard(DM, DM_distr, graphene, "Full", MC_points);
}

}	// namespace Darphene