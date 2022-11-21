#include "graphene/Direct_Detection_NREFT.hpp"

#include <mpi.h>

#include "libphysica/Integration.hpp"
#include "libphysica/List_Manipulations.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"

#include "graphene/Direct_Detection_Standard.hpp"

namespace graphene
{

using namespace libphysica::natural_units;

// 1. Rates and spectra for a single band
// 1.1 Total rate per band
double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, int band, unsigned int MC_points)
{

	// Parameter ranges
	double mDM	= DM.mass;
	double E_b	= graphene.Lowest_Binding_Energy(band);
	double vMin = std::sqrt(2.0 * (graphene.work_function - E_b) / mDM);
	double vMax = DM_distr.Maximum_DM_Speed();
	if(vMax < vMin)
		return 0.0;
	double kfMinGlobal = 0.0;
	double kfMaxGlobal = sqrt(mElectron * mDM) * vMax;
	double qMinGlobal  = mDM * vMax - sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);
	double qMaxGlobal  = mDM * vMax + sqrt(mDM * mDM * vMax * vMax - 2.0 * mDM * graphene.work_function);

	// // NEW INTEGRATION METHOD
	// double prefactor = 1.0 / 32.0 / pow(M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / mDM / mDM / mElectron / mElectron;
	// // Order of integrand arguments: v, phi_v,  xi, cos_theta_q, phi_q, kf, cos_theta_kf, phi_kf
	// std::function<double(const std::vector<double>&, const double)> integrand = [&DM, &DM_distr, &graphene, band, mDM](const std::vector<double>& x, const double wgt) {
	// 	double vDM			= x[0];
	// 	double phi_v		= x[1];
	// 	double xi			= x[2];
	// 	double cos_theta_q	= x[3];
	// 	double phi_q		= x[4];
	// 	double kf			= x[5];
	// 	double cos_theta_kf = x[6];
	// 	double phi_kf		= x[7];

	// 	// 1. Calculate cos_alpha
	// 	// 1.1 kF vector
	// 	Eigen::Vector3d kfVec = Spherical_Coordinates(kf, acos(cos_theta_kf), phi_kf);
	// 	// 1.2 q Vector
	// 	double qMin			 = mDM * vDM - sqrt(mDM * mDM * vDM * vDM - 2.0 * mDM * graphene.work_function);
	// 	double qMax			 = mDM * vDM + sqrt(mDM * mDM * vDM * vDM - 2.0 * mDM * graphene.work_function);
	// 	double q			 = qMin + (qMax - qMin) * xi;
	// 	Eigen::Vector3d qVec = Spherical_Coordinates(q, acos(cos_theta_q), phi_q);

	// 	// 1.3 Calculate kVec ("initial" momentum of the electron)
	// 	Eigen::Vector3d lVec = kfVec - qVec;
	// 	Eigen::Vector3d l_parallel({lVec[0], lVec[1], 0.0});
	// 	Eigen::Vector3d G	 = graphene.Find_G_Vector(l_parallel);
	// 	Eigen::Vector3d kVec = l_parallel - G;

	// 	// 1.4 Calculate band energy
	// 	double E_i = graphene.Valence_Band_Energies(kVec, band);

	// 	// 1.5 Calculate cos_alpha
	// 	double cos_alpha = (kf * kf / 2.0 / mElectron - E_i + q * q / 2.0 / mDM + graphene.work_function) / q / vDM;
	// 	if(cos_alpha < -1.0 || cos_alpha > 1.0)
	// 		return 0.0;

	// 	// 2. vDM vector
	// 	Eigen::Vector3d vDMVec_eigen = Spherical_Coordinates(vDM, acos(cos_alpha), phi_v, qVec);
	// 	libphysica::Vector vDMVec_libphysica({vDMVec_eigen[0], vDMVec_eigen[1], vDMVec_eigen[2]});

	// 	// 3. Evaluate integrand
	// 	return (qMax - qMin) * kf * kf * q * vDM * DM_distr.PDF_Velocity(vDMVec_libphysica) * DM.Response_Function(qVec, vDMVec_eigen, kfVec) * graphene.Material_Response_Function(band, lVec);
	// };
	// double cos_theta_kf_min	   = -1.0;
	// double cos_theta_kf_max	   = 0.0;
	// std::vector<double> region = {vMin, 0.0, 0.0, -1.0, 0.0, kfMinGlobal, cos_theta_kf_min, 0.0, vMax, 2.0 * M_PI, 1.0, 1.0, 2.0 * M_PI, kfMaxGlobal, cos_theta_kf_max, 2.0 * M_PI};
	// double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	// return std::isnan(result) ? 0.0 : result;

	// Order of integrand arguments: q, cos_theta_q, phi_q, cos_theta_kf, phi_kf
	double prefactor = 1.0 / pow(2.0 * M_PI, 2) * DM_distr.DM_density / mDM * graphene.N_cell / 8.0 / mDM / mDM / mElectron / mElectron;

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
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMinGlobal, cos_theta_kf_min, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMaxGlobal, cos_theta_kf_max, 2.0 * M_PI, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 1.2 Differential rates per band
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
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	double cos_theta_kf_min	   = -1.0;
	double cos_theta_kf_max	   = 0.0;
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, cos_theta_kf_min, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, cos_theta_kf_max, 2.0 * M_PI, 1.0, 2.0 * M_PI};
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
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, 0.0, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 2.0 * M_PI, 1.0, 2.0 * M_PI};
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
		Eigen::Vector3d vVec_eigen({v * v_unitvector[0], v * v_unitvector[1], v * v_unitvector[2]});

		return kf * kf * q * v * v / std::fabs(cos_alpha) * DM_distr.PDF_Velocity(vVec) * DM.Response_Function(qVec, vVec_eigen, k_FinalVec) * graphene.Material_Response_Function(band, lVec);
	};
	std::vector<double> region = {qMinGlobal, -1.0, 0.0, kfMin, -1.0, 0.0, qMaxGlobal, 1.0, 2.0 * M_PI, kfMax, 1.0, 2.0 * M_PI};
	double result			   = prefactor * libphysica::Integrate_MC(integrand, region, MC_points, "Vegas");
	return std::isnan(result) ? 0.0 : result;
}

// 2. Rates and spectra for all bands
// 2.1 Total rate
double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double rTot = 0.0;
	for(int band = 0; band < 4; band++)
		rTot += R_Total_NREFT(DM, DM_distr, graphene, band, MC_points);
	return rTot;
}

// 2.2 Differential rates
double dR_dlnE_NREFT(double Ee, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dRdlnE = 0.0;
	for(int band = 0; band < 4; band++)
		dRdlnE += dR_dlnE_NREFT(Ee, DM, DM_distr, graphene, band, MC_points);
	return dRdlnE;
}

double dR_dcos_NREFT(double cos_theta, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_NREFT(cos_theta, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

double dR_dcos_dphi_NREFT(double cos_theta, double phi, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	double dR = 0.0;
	for(int band = 0; band < 4; band++)
		dR += dR_dcos_dphi_NREFT(cos_theta, phi, DM, DM_distr, graphene, band, MC_points);
	return dR;
}

// 3. Tabulation functions
std::vector<std::vector<double>> Tabulate_dR_dlnE_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	// MPI
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	double t_start = MPI_Wtime();

	double E_min			   = 0.1 * eV;
	double v_max			   = DM_distr.Maximum_DM_Speed();
	double E_max			   = 0.99 * DM.mass / 2.0 * v_max * v_max;
	std::vector<double> E_list = libphysica::Log_Space(E_min, E_max, points);
	std::vector<std::vector<double>> local_dRdE(4);
	std::vector<int> index_list = libphysica::Workload_Distribution(mpi_processes, points);

	int counter = 0;
	for(int i = index_list[mpi_rank]; i < index_list[mpi_rank + 1]; i++)
	{
		for(int band = 0; band < 4; band++)
			local_dRdE[band].push_back(dR_dlnE_NREFT(E_list[i], DM, DM_distr, graphene, band, MC_points));
		libphysica::Print_Progress_Bar(1.0 * counter++ / index_list[1], mpi_rank, 50, MPI_Wtime() - t_start);
	}

	// MPI reductions
	std::vector<std::vector<double>> global_dRdE(4, std::vector<double>(points, 0.0));
	std::vector<int> recvcounts(mpi_processes);
	for(int i = 0; i < mpi_processes; i++)
		recvcounts[i] = index_list[i + 1] - index_list[i];
	std::vector<int> displs = libphysica::Sub_List(index_list, 0, mpi_processes);
	for(int band = 0; band < 4; band++)
		MPI_Allgatherv(local_dRdE[band].data(), local_dRdE[band].size(), MPI_DOUBLE, global_dRdE[band].data(), recvcounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	// Compute total dRdE and combine all lists into one.
	std::vector<double> dRdE_total(points, 0.0);
	for(int i = 0; i < points; i++)
		for(int band = 0; band < 4; band++)
			dRdE_total[i] += global_dRdE[band][i];
	std::vector<std::vector<double>> dRdE = libphysica::Transpose_Lists(std::vector<std::vector<double>>({E_list, global_dRdE[0], global_dRdE[1], global_dRdE[2], global_dRdE[3], dRdE_total}));
	libphysica::Print_Progress_Bar(1.0, mpi_rank, 50, MPI_Wtime() - t_start);
	return dRdE;
}

std::vector<std::vector<double>> Tabulate_dR_dcos_dphi_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	// MPI
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	double t_start = MPI_Wtime();

	auto cos_k_list = libphysica::Linear_Space(-1.0, 1.0, points);
	auto phi_k_list = libphysica::Linear_Space(0.0, 2.0 * M_PI, points);

	std::vector<std::vector<double>> local_output(3);
	std::vector<int> index_list = libphysica::Workload_Distribution(mpi_processes, points);
	int counter					= 0;

	for(int i = index_list[mpi_rank]; i < index_list[mpi_rank + 1]; i++)
	{
		double cos_theta = cos_k_list[i];
		for(auto& phi : phi_k_list)
		{
			local_output[0].push_back(cos_theta);
			local_output[1].push_back(phi);
			local_output[2].push_back(dR_dcos_dphi_NREFT(cos_theta, phi, DM, DM_distr, graphene, MC_points));
			libphysica::Print_Progress_Bar(1.0 * counter++ / index_list[1] / points, mpi_rank, 50, MPI_Wtime() - t_start);
		}
	}

	// Gather all rates from the MPI processes.
	std::vector<std::vector<double>> global_output(3, std::vector<double>(points * points, 0.0));
	std::vector<int> recvcounts(mpi_processes);
	for(int i = 0; i < mpi_processes; i++)
		recvcounts[i] = (index_list[i + 1] - index_list[i]) * points;
	std::vector<int> displs = libphysica::Sub_List(index_list, 0, mpi_processes);
	for(auto& displ : displs)
		displ *= points;
	for(int i = 0; i < 3; i++)
		MPI_Allgatherv(local_output[i].data(), local_output[i].size(), MPI_DOUBLE, global_output[i].data(), recvcounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	libphysica::Print_Progress_Bar(1.0, mpi_rank, 50, MPI_Wtime() - t_start);
	return libphysica::Transpose_Lists(global_output);
}

std::vector<std::vector<double>> Daily_Modulation_NREFT(int points, DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points)
{
	// MPI
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	double t_start = MPI_Wtime();

	// Save original observer velocity
	libphysica::Vector observer_velocity = dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Get_Observer_Velocity();
	double vEarth						 = observer_velocity.Norm();

	// Total rate over the course of a day
	std::vector<double> t_list = libphysica::Linear_Space(0.0, 24.0, points);
	std::vector<double> local_rates;
	std::vector<int> index_list = libphysica::Workload_Distribution(mpi_processes, points);
	int counter					= 0;
	for(int i = index_list[mpi_rank]; i < index_list[mpi_rank + 1]; i++)
	{
		// Set observer velocity at given time t.
		double t = t_list[i];
		dynamic_cast<obscura::Standard_Halo_Model*>(&DM_distr)->Set_Observer_Velocity(Earth_Velocity(t * hr, vEarth));

		double R = R_Total_NREFT(DM, DM_distr, graphene, MC_points);
		local_rates.push_back(R);
		libphysica::Print_Progress_Bar(1.0 * counter++ / index_list[1], mpi_rank, 50, MPI_Wtime() - t_start);
	}

	// Gather all rates from the MPI processes.
	std::vector<double> global_rates(points);
	std::vector<int> recvcounts(mpi_processes);
	for(int i = 0; i < mpi_processes; i++)
		recvcounts[i] = index_list[i + 1] - index_list[i];
	std::vector<int> displs = libphysica::Sub_List(index_list, 0, mpi_processes);

	MPI_Allgatherv(local_rates.data(), local_rates.size(), MPI_DOUBLE, global_rates.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	std::vector<std::vector<double>> daily_modulation_list = libphysica::Transpose_Lists(t_list, global_rates);

	libphysica::Print_Progress_Bar(1.0, mpi_rank, 50, MPI_Wtime() - t_start);
	return daily_modulation_list;
}

}	// namespace graphene