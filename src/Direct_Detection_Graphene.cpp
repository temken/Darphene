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
	return (E_final + energy_crystal + work_function) / q + q / 2.0 / mDM;
}

double dR_dlnE(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, const std::string& method)
{
	//1. Prefactor
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
	std::function<double(const std::vector<double>&, const double)> integrand = [&DM_distr, &graphene, &DM, band, k_final, work_function, E_e, mDM, vMax, prefactor](const std::vector<double>& x, const double wgt) {
		double l_1			= x[0];
		double xi			= x[1];
		double y			= x[2];
		double cos_theta_q	= x[3];
		double phi_q		= x[4];
		double cos_theta_kf = x[5];
		double phi_kf		= x[6];

		double l_2 = l_1 / sqrt(3.0) * xi;
		Eigen::Vector3d lVec({l_1, l_2, 0.0});

		double E_l				   = graphene.Valence_Band_Energies(lVec, band);
		double qMin				   = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_l + work_function + E_e));
		double qMax				   = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_l + work_function + E_e));
		double q				   = qMin + y * (qMax - qMin);
		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(k_final, acos(cos_theta_kf), phi_kf);
		double vMin				   = vMinimum_Graphene(DM.mass, q, E_l, k_final);
		double vDM				   = 1.0e-3;   // cancels
		return lVec[0] * (qMax - qMin) * q * DM_distr.Eta_Function(vMin) * vDM * vDM * DM.dSigma_dq2_Electron(q, vDM) * graphene.DM_Response(band, lVec, qVec - k_FinalVec);
	};
	std::vector<double> region = {0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 2.0 * M_PI / sqrt(3.0) / a, 1.0, 1.0, +1.0, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	return prefactor * libphysica::Integrate_MC(integrand, region, 10000, method);
}

}	// namespace graphene