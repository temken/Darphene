#include "Direct_Detection_Graphene.hpp"

#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

#include "Multidimensional_Integration.hpp"
#include "Vegas.hpp"

namespace graphene
{

using namespace libphysica::natural_units;

double vMinimum_Graphene(double mDM, double q, double energy_crystal, double final_momentum)
{
	double E_final		 = final_momentum * final_momentum / 2.0 / mElectron;
	double work_function = 4.3 * eV;
	return (E_final + energy_crystal + work_function) / q + q / 2.0 / mDM;
}

double perform_integral(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band)
{
	std::random_device rd;
	std::mt19937 PRNG(rd());

	//1. Prefactor
	double mDM		 = DM.mass;
	double mu		 = libphysica::Reduced_Mass(mDM, mElectron);
	double sigma_e	 = 0.1 * pb;
	double aCC		 = 1.42 * Angstrom;
	double a		 = aCC * sqrt(3.0);
	double N_C		 = 5.0e25 / kg;
	double A_UC		 = 3.0 * sqrt(3.0) * aCC * aCC / 2.0;
	double k_final	 = sqrt(2.0 * mElectron * E_e);
	double prefactor = sigma_e * DM_distr.DM_density / mDM * N_C * A_UC * mElectron * k_final / pow(2.0 * M_PI, 5.0) / mu / mu * E_e;

	double work_function = 4.3 * eV;
	double vMax			 = DM_distr.Maximum_DM_Speed();

	unsigned int Nmax = 100000;
	double sum		  = 0.0;
	double volume	  = 128.0 * M_PI * M_PI * M_PI * M_PI / sqrt(3.0) / a / a;
	for(int i = 0; i < Nmax; i++)
	{
		// 1. Sample l in the (G,M,K) triangle of the 1st BZ
		double xi_1 = libphysica::Sample_Uniform(PRNG);
		double xi_2 = libphysica::Sample_Uniform(PRNG);
		if(xi_1 + xi_2 > 1)
		{
			xi_1 = 1.0 - xi_1;
			xi_2 = 1.0 - xi_2;
		}
		Eigen::Vector3d lVec({(xi_1 + xi_2) * 2.0 * M_PI / sqrt(3.0) / a, xi_2 * 2.0 * M_PI / 3.0 / a, 0.0});

		// 2. Sample q
		double E_crystal	 = graphene.Valence_Band_Energies(lVec, band);
		double qMin			 = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_crystal + work_function + E_e));
		double qMax			 = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_crystal + work_function + E_e));
		double q			 = libphysica::Sample_Uniform(PRNG, qMin, qMax);
		double cos_theta_q	 = libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		double phi_q		 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
		Eigen::Vector3d qVec = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);

		//3. Sample direction of k_final
		double cos_theta_k		  = libphysica::Sample_Uniform(PRNG, -1.0, 1.0);
		double phi_k			  = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
		Eigen::Vector3d kFinalVec = Spherical_Coordinates(k_final, acos(cos_theta_k), phi_k);
		// sum += (qMax - qMin) * integrand(band, DM, DM_distr, graphene, E_crystal, lVec, qVec, kFinalVec);
		double F_DM = 1.0;
		sum += (qMax - qMin) * q * DM_distr.Eta_Function(vMinimum_Graphene(DM.mass, qVec.norm(), E_crystal, k_final)) * F_DM * graphene.DM_Response(band, lVec, qVec - kFinalVec);
	}
	double integral = volume / Nmax * sum;
	return prefactor * integral;
}

double perform_integral_vegas(double E_e, obscura::DM_Particle& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band)
{
	//1. Prefactor
	double mDM		 = DM.mass;
	double mu		 = libphysica::Reduced_Mass(mDM, mElectron);
	double sigma_e	 = 0.1 * pb;
	double aCC		 = 1.42 * Angstrom;
	double a		 = aCC * sqrt(3.0);
	double N_C		 = 5.0e25 / kg;
	double A_UC		 = 3.0 * sqrt(3.0) * aCC * aCC / 2.0;
	double k_final	 = sqrt(2.0 * mElectron * E_e);
	double prefactor = 12.0 * sigma_e * DM_distr.DM_density / mDM * N_C * A_UC * mElectron * k_final / pow(2.0 * M_PI, 5.0) / mu / mu * E_e;

	double work_function = 4.3 * eV;
	double vMax			 = DM_distr.Maximum_DM_Speed();

	// 2. Integrate with VEGAS algorithm
	// Order of integrand arguments: l_1, xi = l_2 / (l_1/sqrt(3)) , y = (q-qMin)/(qMax-qMin), cos_theta_q, phi_q, cos_theta_kf, phi_kf
	std::vector<double> regn = {0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 2.0 * M_PI / sqrt(3.0) / a, 1.0, 1.0, +1.0, 2.0 * M_PI, 0.0, 2.0 * M_PI};

	std::function<double(const std::vector<double>&, const double)> integrand_vegas = [&DM_distr, &graphene, &DM, band, k_final, work_function, E_e, mDM, vMax](const std::vector<double>& x, const double wgt) {
		double l_1 = x[0];
		double xi  = x[1];
		double l_2 = l_1 / sqrt(3.0) * xi;
		Eigen::Vector3d lVec({l_1, l_2, 0.0});
		double y			= x[2];
		double cos_theta_q	= x[3];
		double phi_q		= x[4];
		double cos_theta_kf = x[5];
		double phi_kf		= x[6];

		double E_l				   = graphene.Valence_Band_Energies(lVec, band);
		double qMin				   = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_l + work_function + E_e));
		double qMax				   = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * (E_l + work_function + E_e));
		double q				   = qMin + y * (qMax - qMin);
		double F_DM				   = 1.0;
		Eigen::Vector3d qVec	   = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);
		Eigen::Vector3d k_FinalVec = Spherical_Coordinates(k_final, acos(cos_theta_kf), phi_kf);
		double vMin				   = vMinimum_Graphene(DM.mass, q, E_l, k_final);
		return lVec[0] / sqrt(3.0) * (qMax - qMin) * q * DM_distr.Eta_Function(vMin) * F_DM * graphene.DM_Response(band, lVec, qVec - k_FinalVec);
	};
	double tgral, sd, chi2a;
	vegas(regn, integrand_vegas, 0, 100000, 10, -1, tgral, sd, chi2a);
	return prefactor * tgral;
}

}	// namespace graphene