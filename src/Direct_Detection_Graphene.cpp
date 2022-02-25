#include "Direct_Detection_Graphene.hpp"

#include "obscura/DM_Halo_Models.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

namespace graphene
{

using namespace libphysica::natural_units;

Eigen::Vector3d Spherical_Coordinates(double r, double theta, double phi)
{
	return Eigen::Vector3d(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
}

extern libphysica::Vector Earth_Velocity(double t, double v_earth)
{
	double alpha = 42.0 * deg;
	double beta	 = 2.0 * M_PI * t / day;

	double cosa = cos(alpha);
	double sina = sin(alpha);
	double cosb = cos(beta);
	double sinb = sin(beta);

	return libphysica::Vector({v_earth * sina * sinb, v_earth * sina * cosa * (cosb - 1.0), v_earth * (cosa * cosa + sina * sina * cosb)});
}

Eigen::Vector3d Spherical_Coordinates_Axis(double r, double theta, double phi, const libphysica::Vector& axis)
{
	libphysica::Vector ev = axis.Normalized();
	if(ev[2] == 1.0 || axis.Norm() == 0.0)
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

double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum, double work_function)
{
	double E_final = final_momentum * final_momentum / 2.0 / mElectron;
	return (E_final - energy_crystal + work_function) / q + q / 2.0 / mDM;
}

double dR_dlnE(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double aCC		 = 1.42 * Angstrom;
	double a		 = aCC * sqrt(3.0);
	double N_C		 = 5.0e25 / kg;
	double A_UC		 = 3.0 * sqrt(3.0) * aCC * aCC / 2.0;
	double k_final	 = sqrt(2.0 * mElectron * E_e);
	double prefactor = 8.0 * sqrt(3.0) * DM_distr.DM_density / mDM * N_C * A_UC * mElectron * k_final / pow(2.0 * M_PI, 6.0) * E_e;

	double vMax = DM_distr.Maximum_DM_Speed();

	// Order of integrand arguments: l_1, xi = l_2 / (l_1/sqrt(3)) , y = (q-qMin)/(qMax-qMin), cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, k_final, E_e, mDM, vMax](const std::vector<double>& x, const double wgt) {
		double l_1			= x[0];
		double xi			= x[1];
		double y			= x[2];
		double cos_theta_q	= x[3];
		double phi_q		= x[4];
		double cos_theta_kf = x[5];
		double phi_kf		= x[6];

		double l_2 = l_1 / sqrt(3.0) * xi;
		Eigen::Vector3d lVec({l_1, l_2, 0.0});

		double E_l				   = graphene.Valence_Band_Energies(lVec, band);   // Note that E_l is negative!!! (unlike in the paper by Hochberg et al)
		double qMin				   = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function - E_l + E_e));
		double qMax				   = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function - E_l + E_e));
		double q				   = qMin + y * (qMax - qMin);
		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(k_final, acos(cos_theta_kf), phi_kf);
		double vMin				   = vMinimum_Graphene(DM.mass, q, E_l, k_final, graphene.work_function);
		double vDM				   = 1.0e-3;   // cancels
		return lVec[0] * (qMax - qMin) * q * DM_distr.Eta_Function(vMin) * vDM * vDM * DM.dSigma_dq2_Electron(q, vDM) * graphene.DM_Response(band, lVec, qVec - k_FinalVec);
	};
	std::vector<double> region = {0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 2.0 * M_PI / sqrt(3.0) / a, 1.0, 1.0, +1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 5000, method);
	return std::isnan(result) ? 0.0 : result;
}

double dR_dlnE_corrected(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double N_C		 = 5.0e25 / kg;
	double k_final	 = sqrt(2.0 * mElectron * E_e);
	double prefactor = 0.5 / pow(2.0 * M_PI, 4) * DM_distr.DM_density / mDM * N_C * sqrt(2.0 * pow(mElectron * E_e, 3.0)) * sigma_e / mu_e / mu_e;

	double vMax = DM_distr.Maximum_DM_Speed();
	double EMax = mDM / 2.0 * vMax * vMax;
	if(E_e > EMax)
		return 0.0;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function + E_e));
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function + E_e));

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, k_final, E_e, mDM, vMax, qMinGlobal, qMaxGlobal](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double cos_theta_kf = x[3];
		double phi_kf		= x[4];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(k_final, acos(cos_theta_kf), phi_kf);

		// Determine the crystal momentum vector l
		Eigen::Vector3d lVec({k_FinalVec[0] - qVec[0], k_FinalVec[1] - qVec[1], 0.0});
		lVec = graphene.Find_1BZ_Vector(lVec);

		double E_l = graphene.Valence_Band_Energies(lVec, band);   // Note that E_l is negative!!! (unlike in the paper by Hochberg et al)

		double vMin = vMinimum_Graphene(DM.mass, q, E_l, k_final, graphene.work_function);
		return q * DM_distr.Eta_Function(vMin) * graphene.DM_Response_corrected(band, qVec, k_FinalVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 5000, method);

	return std::isnan(result) ? 0.0 : result;
}

double R_Total_corrected(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double N_C		 = 5.0e25 / kg;
	double prefactor = 0.5 / pow(2.0 * M_PI, 4) * DM_distr.DM_density / mDM * N_C * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, mDM, vMax, qMinGlobal, qMaxGlobal, kfMin, kfMax](const std::vector<double>& x, const double wgt) {
		double q			= x[0];
		double cos_theta_q	= x[1];
		double phi_q		= x[2];
		double kf			= x[3];
		double cos_theta_kf = x[4];
		double phi_kf		= x[5];

		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);

		// Determine the crystal momentum vector l
		Eigen::Vector3d lVec({k_FinalVec[0] - qVec[0], k_FinalVec[1] - qVec[1], 0.0});
		lVec = graphene.Find_1BZ_Vector(lVec);

		double E_l = graphene.Valence_Band_Energies(lVec, band);   // Note that E_l is negative!!! (unlike in the paper by Hochberg et al)

		double vMin = vMinimum_Graphene(DM.mass, q, E_l, kf, graphene.work_function);
		return kf * kf * q * DM_distr.Eta_Function(vMin) * graphene.DM_Response_corrected(band, qVec, k_FinalVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 200000, method);

	return std::isnan(result) ? 0.0 : result;
}
double R_Total_corrected(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method)
{
	double rTot = 0.0;
	for(int band = 0; band < 4; band++)
		rTot += R_Total_corrected(DM, DM_distr, graphene, band, method);
	return rTot;
}

double R_Total_Full_Integral(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double N_C		 = 5.0e25 / kg;
	double prefactor = 1.0 / pow(2.0 * M_PI, 4) * DM_distr.DM_density / mDM * N_C * sigma_e / mu_e / mu_e;

	double vMax	 = DM_distr.Maximum_DM_Speed();
	double kfMin = 0.0;
	double kfMax = sqrt(mElectron * mDM) * vMax;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector l
		Eigen::Vector3d lVec({k_FinalVec[0] - qVec[0], k_FinalVec[1] - qVec[1], 0.0});
		lVec = graphene.Find_1BZ_Vector(lVec);

		double E_l = graphene.Valence_Band_Energies(lVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_l + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.DM_Response_corrected(band, qVec, k_FinalVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 100000, method);

	return std::isnan(result) ? 0.0 : result;
}

double R_Total_Full_Integral(obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, const std::string& method)
{
	double rTot = 0.0;
	for(int band = 0; band < 4; band++)
		rTot += R_Total_Full_Integral(DM, DM_distr, graphene, band, method);
	return rTot;
}

double dR_dlnE_Full_Integral(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	// 1. Prefactor
	double mDM		 = DM.mass;
	double sigma_e	 = DM.Sigma_Electron();
	double mu_e		 = libphysica::Reduced_Mass(mElectron, mDM);
	double N_C		 = 5.0e25 / kg;
	double prefactor = std::sqrt(2.0 * std::pow(mElectron * E_e, 3.0)) / pow(2.0 * M_PI, 4) * DM_distr.DM_density / mDM * N_C * sigma_e / mu_e / mu_e;

	double kf = std::sqrt(2.0 * mElectron * E_e);

	double vMax = DM_distr.Maximum_DM_Speed();
	double EMax = mDM / 2.0 * vMax * vMax;
	if(E_e > EMax)
		return 0.0;

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function + E_e));
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (graphene.work_function + E_e));

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [kf, &DM_distr, &graphene, &DM, band, mDM, vMax](const std::vector<double>& x, const double wgt) {
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

		// Determine the crystal momentum vector l
		Eigen::Vector3d lVec({k_FinalVec[0] - qVec[0], k_FinalVec[1] - qVec[1], 0.0});
		lVec = graphene.Find_1BZ_Vector(lVec);

		double E_l = graphene.Valence_Band_Energies(lVec, band);
		double v   = (kf * kf / (2.0 * mElectron) - E_l + graphene.work_function + q * q / 2.0 / mDM) / (q * cos_alpha);
		if(v > vMax || v < 0.0)
			return 0.0;
		libphysica::Vector vVec({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * graphene.DM_Response_corrected(band, qVec, k_FinalVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 20000, method);

	return std::isnan(result) ? 0.0 : result;
}

}	// namespace graphene