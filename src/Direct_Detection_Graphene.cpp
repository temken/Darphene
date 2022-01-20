#include "Direct_Detection_Graphene.hpp"

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

double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum)
{
	double E_final		 = final_momentum * final_momentum / 2.0 / mElectron;
	double work_function = 4.3 * eV;
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

	double work_function = 4.3 * eV;
	double vMax			 = DM_distr.Maximum_DM_Speed();

	// Order of integrand arguments: l_1, xi = l_2 / (l_1/sqrt(3)) , y = (q-qMin)/(qMax-qMin), cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, k_final, work_function, E_e, mDM, vMax](const std::vector<double>& x, const double wgt) {
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
		double qMin				   = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (work_function - E_l + E_e));
		double qMax				   = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (work_function - E_l + E_e));
		double q				   = qMin + y * (qMax - qMin);
		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(k_final, acos(cos_theta_kf), phi_kf);
		double vMin				   = vMinimum_Graphene(DM.mass, q, E_l, k_final);
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

	double work_function = 4.3 * eV;
	double vMax			 = DM_distr.Maximum_DM_Speed();

	double qMinGlobal = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (work_function + E_e));
	double qMaxGlobal = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (work_function + E_e));

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, k_final, work_function, E_e, mDM, vMax, qMinGlobal, qMaxGlobal](const std::vector<double>& x, const double wgt) {
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

		double vMin = vMinimum_Graphene(DM.mass, q, E_l, k_final);
		return q * DM_distr.Eta_Function(vMin) * graphene.DM_Response_corrected(band, qVec, k_FinalVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, 5000, method);

	return std::isnan(result) ? 0.0 : result;
}

}	// namespace graphene