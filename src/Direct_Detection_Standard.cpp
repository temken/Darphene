#include "graphene/Direct_Detection_Standard.hpp"

#include <omp.h>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"

namespace graphene
{

using namespace libphysica::natural_units;

Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi)
{
	return Eigen::Vector3d(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
}

// 1. Total scattering rates
// 1.1 Total rate per band
double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");

	return std::isnan(result) ? 0.0 : result;
}

double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum, double work_function)
{
	double E_final = final_momentum * final_momentum / 2.0 / mElectron;
	return (E_final - energy_crystal + work_function) / q + q / 2.0 / mDM;
}

double R_Total_Standard_Simple(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k	= graphene.Valence_Band_Energies(kVec, band);
		double vMin = vMinimum_Graphene(mDM, q, E_k, kf, graphene.work_function);

		return kf * kf * q * DM_distr.Eta_Function(vMin) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI};

	double result = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double prefactor = 1.0 / pow(2.0 * M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / 8.0 / mDM / mDM / mElectron / mElectron;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM, &DM_distr, &graphene, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		// double E_k = graphene.Valence_Band_Energies(kVec, band);
		double E_k = -9.0 * eV;
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 1.2 Total rate for all bands
double R_Total_Standard(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points)
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
double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double rTot = 0.0;
	for(int band = 0; band < 4; band++)
		rTot += R_Total_NREFT(DM, DM_distr, graphene, band, MC_points);
	return rTot;
}

// 2. Differential rate dR/dlnE
// 2.1 For one band
double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

double dR_dlnE_Standard_Simple(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k	= graphene.Valence_Band_Energies(kVec, band);
		double vMin = vMinimum_Graphene(mDM, q, E_k, kf, graphene.work_function);

		return q * DM_distr.Eta_Function(vMin) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

double dR_dlnE_NREFT(double Ee, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double prefactor = std::sqrt(2.0 * std::pow(mElectron * Ee, 3.0)) / pow(2.0 * M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / 8.0 / mDM / mDM / mElectron / mElectron;

	double kf = std::sqrt(2.0 * mElectron * Ee);

	double vMax = DM_distr.Maximum_DM_Speed();
	double EMax = mDM / 2.0 * vMax * vMax;
	if(Ee + graphene.work_function > EMax)
		return 0.0;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM, &DM_distr, &graphene, band, mDM, vMax, kf](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 2.2 For all bands
double dR_dlnE_Standard(double Ee, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points)
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

double dR_dlnE_NREFT(double Ee, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dRdlnE = 0.0;
	for(int band = 0; band < 4; band++)
		dRdlnE += dR_dlnE_NREFT(Ee, DM, DM_distr, graphene, band, MC_points);
	return dRdlnE;
}

// 3. Differential rate dR / dcos(theta)
// 3.1 For one band
double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

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
double dR_dcos_NREFT(double cos_theta, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double prefactor = 1.0 / pow(2.0 * M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / 8.0 / mDM / mDM / mElectron / mElectron;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [cos_theta, &DM, &DM_distr, &graphene, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		// double E_k = graphene.Valence_Band_Energies(kVec, band);
		double E_k = -9.0 * eV;
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 3.2 For all bands
double dR_dcos_Standard(double cos_theta, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_Standard(cos_theta, DM, DM_distr, graphene, band, MC_points);
	return dR;
}
double dR_dcos_NREFT(double cos_theta, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_NREFT(cos_theta, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

// 4. Differential rate d^2R/(dcos dphi)
// 4.1 For one band
double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

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

double dR_dcos_dphi_NREFT(double cos_theta, double phi, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double prefactor = 1.0 / pow(2.0 * M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / 8.0 / mDM / mDM / mElectron / mElectron;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM, &DM_distr, &graphene, band, mDM, vMax, cos_theta, phi](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector k = G* - l^|| (in 1BZ)
		Eigen::Vector3d lVec = k_FinalVec - qVec;
		Eigen::Vector3d kVec({-lVec[0], -lVec[1], 0.0});
		kVec = graphene.Find_1BZ_Vector(kVec);

		double E_k = graphene.Valence_Band_Energies(kVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_k + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 4.2 For all bands
double dR_dcos_dphi_Standard(double cos_theta, double phi, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_dphi_Standard(cos_theta, phi, DM, DM_distr, graphene, band, MC_points);
	return dR;
}
double dR_dcos_dphi_NREFT(double cos_theta, double phi, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_dphi_NREFT(cos_theta, phi, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

// 5. Tabulation functions
// 5.1 Energy spectrum
std::vector<std::vector<double>> Tabulate_dR_dlnE_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& velocity_integral, unsigned int MC_points, int threads)
{
	double E_min			   = 0.1 * eV;
	double v_max			   = DM_distr.Maximum_DM_Speed();
	double E_max			   = 0.99 * DM.mass / 2.0 * v_max * v_max;
	std::vector<double> E_kist = libphysica::Log_Space(E_min, E_max, points);
	std::vector<std::vector<double>> spectrum;
#pragma omp parallel for schedule(dynamic) num_threads(threads)
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

std::vector<std::vector<double>> Tabulate_dR_dlnE_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads)
{
	double E_min			   = 0.1 * eV;
	double v_max			   = DM_distr.Maximum_DM_Speed();
	double E_max			   = 0.99 * DM.mass / 2.0 * v_max * v_max;
	std::vector<double> E_kist = libphysica::Log_Space(E_min, E_max, points);
	std::vector<std::vector<double>> spectrum;
	int counter = 0;
#pragma omp parallel for schedule(dynamic) num_threads(threads)
	for(auto& E_er : E_kist)
	{
		std::vector<double> row = {E_er};
		double dRdlnE			= 0.0;
		for(int band = 0; band < 4; band++)
		{
			double band_contribution = dR_dlnE_NREFT(E_er, DM, DM_distr, graphene, band, MC_points);
			row.push_back(band_contribution);
			dRdlnE += band_contribution;
		}
		row.push_back(dRdlnE);
		spectrum.push_back(row);
		libphysica::Print_Progress_Bar(1.0 * counter++ / points);
	}
	libphysica::Print_Progress_Bar(1.0);

	return spectrum;
}

// 5.2 Directional spectrum
std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads)
{
	auto cos_k_list = libphysica::Linear_Space(-1.0, 1.0, points);
	auto phi_k_list = libphysica::Linear_Space(0.0, 2.0 * M_PI, points);

	std::vector<std::vector<double>> spectrum;
#pragma omp parallel for schedule(dynamic) num_threads(threads) collapse(2)
	for(auto& cos_theta : cos_k_list)
		for(auto& phi : phi_k_list)
			spectrum.push_back({cos_theta, phi, dR_dcos_dphi_Standard(cos_theta, phi, DM, DM_distr, graphene, MC_points)});

	return spectrum;
}

std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads)
{
	auto cos_k_list = libphysica::Linear_Space(-1.0, 1.0, points);
	auto phi_k_list = libphysica::Linear_Space(0.0, 2.0 * M_PI, points);

	std::vector<std::vector<double>> spectrum;
	int counter = 0;
#pragma omp parallel for schedule(dynamic) num_threads(threads) collapse(2)
	for(auto& cos_theta : cos_k_list)
		for(auto& phi : phi_k_list)
		{
			spectrum.push_back({cos_theta, phi, dR_dcos_dphi_NREFT(cos_theta, phi, DM, DM_distr, graphene, MC_points)});
			libphysica::Print_Progress_Bar(1.0 * counter++ / points / points);
		}
	libphysica::Print_Progress_Bar(1.0);
	return spectrum;
}

// 5.3 Daily Modulation
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

std::vector<std::vector<double>> Daily_Modulation_Standard(int points, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads)
{
	// Total rate over the course of a day
	std::vector<double> t_list = libphysica::Linear_Space(0.0, 24.0, points);
	std::vector<std::vector<double>> daily_modulation_list;
	double vEarth = dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Get_Observer_Velocity().Norm();
#pragma omp parallel for schedule(dynamic) num_threads(threads)
	for(auto& t : t_list)
	{
		dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Set_Observer_Velocity(Earth_Velocity(t * hr, vEarth));
		double R = R_Total_Standard(DM, DM_distr, graphene, "Full", MC_points);
		daily_modulation_list.push_back({t, R});
	}
	return daily_modulation_list;
}

std::vector<std::vector<double>> Daily_Modulation_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points, int threads)
{
	// Total rate over the course of a day
	std::vector<double> t_list = libphysica::Linear_Space(0.0, 24.0, points);
	std::vector<std::vector<double>> daily_modulation_list;
	double vEarth = dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Get_Observer_Velocity().Norm();
	int counter	  = 0;
#pragma omp parallel for schedule(dynamic) num_threads(threads)
	for(auto& t : t_list)
	{
		libphysica::Print_Progress_Bar(1.0 * counter++ / points);
		dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Set_Observer_Velocity(Earth_Velocity(t * hr, vEarth));
		double R = R_Total_NREFT(DM, DM_distr, graphene, MC_points);
		daily_modulation_list.push_back({t, R});
	}
	libphysica::Print_Progress_Bar(1.0);
	return daily_modulation_list;
}

}	// namespace graphene