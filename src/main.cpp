#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>
#include <random>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

#include "graphene/Carbon_Wavefunctions.hpp"
#include "graphene/Configuration.hpp"
#include "graphene/DM_Particle_NREFT.hpp"
#include "graphene/Direct_Detection_NREFT.hpp"
#include "graphene/Direct_Detection_Standard.hpp"
#include "graphene/Graphene.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace graphene;

int main(int argc, char* argv[])
{
	// MPI initialization
	MPI_Init(NULL, NULL);
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	if(mpi_rank == 0)
	{
		std::cout << "[Started on " << ctime_start << "]" << std::endl;
		std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
				  << std::endl
				  << "MPI processes:\t" << mpi_processes << std::endl;
	}
	////////////////////////////////////////////////////////////////////////

	Configuration cfg(argv[1]);
	Graphene graphene(cfg.carbon_wavefunctions, cfg.graphene_work_function);
	cfg.Print_Summary(mpi_rank);
	std::string results_path = TOP_LEVEL_DIR "results/" + cfg.ID + "/";
	double rate_unit		 = 1.0 / gram / year;

	MPI_Barrier(MPI_COMM_WORLD);
	if(cfg.run_modus == "Energy-Spectrum" || cfg.run_modus == "All")
	{
		if(mpi_rank == 0)
			std::cout << "\nTabulate dR/dlnE:" << std::endl;
		auto spectrum_nreft = Tabulate_dR_dlnE_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "dR_dlnE.txt";
			libphysica::Export_Table(file_path, spectrum_nreft, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});
			std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
		}
	}
	if(cfg.run_modus == "Total-Rate" || cfg.run_modus == "All")
	{
		if(mpi_rank == 0)
			std::cout << "\nCompute total rate R:" << std::endl;
		double R = R_Total_NREFT(*cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Total_Rate.txt";
			std::ofstream f(file_path);
			f << "R = " << In_Units(R, rate_unit) << " /gr / year" << std::endl;
			std::cout << "R = " << In_Units(R, rate_unit) << " / gr / year" << std::endl;
			f.close();
			std::cout << "\nDone. Total rate saved to " << file_path << "." << std::endl;
		}
	}
	if(cfg.run_modus == "Directional-Spectrum" || cfg.run_modus == "All")
	{
		if(mpi_rank == 0)
			std::cout << "\nTabulate dR/(dcos dphi):" << std::endl;
		auto spectrum_nreft = Tabulate_dR_dcos_dphi_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "dR_dcos_dphi.txt";
			libphysica::Export_Table(file_path, spectrum_nreft, {1.0, 1.0, rate_unit}, "# cos_theta\tphi [rad]\tdR/(dcos dphi) [1/yr/gr]");
			std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
		}
	}
	if(cfg.run_modus == "Daily-Modulation" || cfg.run_modus == "All")
	{
		if(mpi_rank == 0)
			std::cout << "\nCalculate daily modulation:" << std::endl;
		auto daily_nreft = Daily_Modulation_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Daily_Modulation.txt";
			libphysica::Export_Table(file_path, daily_nreft, {1.0, rate_unit}, "# t[h]\tR[1/gr/year]");
			std::cout << "\nDone. Table saved to " << file_path << "." << std::endl;
		}
	}
	else if(cfg.run_modus == "Custom")
	{
		// for(auto overlap : {"s", "sPrime", "Sss", "Ssigma", "Ssp"})
		// 	std::cout << overlap << "\t" << graphene.Overlap_Integral(overlap) << std::endl;

		// // 1. DM model and mass
		// // std::vector<double> DM_masses = libphysica::Log_Space(mMin, 100.0 * MeV, 10);
		// std::vector<int> DM_masses = {2, 5, 10, 20, 50, 100};

		// std::vector<int> operators			  = {1, 3};
		// std::vector<std::string> interactions = {"C", "L"};
		// for(auto& op : operators)
		// 	for(auto& inter : interactions)
		// 		for(auto& mDM : DM_masses)
		// 		{
		// 			if(mpi_rank == 0)
		// 				std::cout << "Operator " << op << "\tInteraction " << inter << "\tDM mass = " << mDM << " MeV" << std::endl;
		// 			DM_Particle_NREFT DM(mDM * MeV);
		// 			if(inter == "C")
		// 				DM.Set_Coupling(op, 1.0, "Contact");
		// 			else if(inter == "L")
		// 				DM.Set_Coupling(op, 1.0, "Long-Range");
		// 			// 2. Compute the daily modulation
		// 			std::string file_path = results_path + "Daily_Modulation_O" + std::to_string(op) + inter + "_m_" + std::to_string(mDM) + "MeV.txt";
		// 			if(libphysica::File_Exists(file_path) == false)
		// 			{
		// 				std::vector<std::vector<double>> daily_modulation = Daily_Modulation_NREFT(25, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 				if(mpi_rank == 0)
		// 					libphysica::Export_Table(file_path, daily_modulation, {1.0, 1.0 / gram / year}, "#Time [h]\tRate [1/gr/year]");
		// 			}
		// 		}

		// for(auto& mDM : DM_masses)
		// {
		// 	std::cout << "Electric dipole with mDM = " << mDM << " MeV" << std::endl;
		// 	DM_Particle_NREFT DM = DM_Electric_Dipole(mDM * MeV, 1.0);
		// 	// 2. Compute the daily modulation
		// 	std::string file_path = results_path + "Daily_Modulation_ED_m_" + std::to_string(mDM) + "MeV.txt";
		// 	if(libphysica::File_Exists(file_path) == false)
		// 	{
		// 		std::vector<std::vector<double>> daily_modulation = Daily_Modulation_NREFT(25, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 		if(mpi_rank == 0)
		// 			libphysica::Export_Table(file_path, daily_modulation, {1.0, 1.0 / gram / year}, "#Time [h]\tRate [1/gr/year]");
		// 	}
		// }

		// for(auto& mDM : DM_masses)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "Magnetic dipole with mDM = " << mDM << " MeV" << std::endl;
		// 	DM_Particle_NREFT DM = DM_Magnetic_Dipole(mDM * MeV, 1.0);

		// 	// 2. Compute the daily modulation
		// 	std::string file_path = results_path + "Daily_Modulation_MD_m_" + std::to_string(mDM) + "MeV.txt";
		// 	if(libphysica::File_Exists(file_path) == false)
		// 	{
		// 		std::vector<std::vector<double>> daily_modulation = Daily_Modulation_NREFT(25, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 		if(mpi_rank == 0)
		// 			libphysica::Export_Table(file_path, daily_modulation, {1.0, 1.0 / gram / year}, "#Time [h]\tRate [1/gr/year]");
		// 	}
		// }

		// for(auto& mDM : DM_masses)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "Anapole with mDM = " << mDM << " MeV" << std::endl;
		// 	DM_Particle_NREFT DM = DM_Anapole(mDM * MeV, 1.0);

		// 	// 2. Compute the daily modulation
		// 	std::string file_path = results_path + "Daily_Modulation_A_m_" + std::to_string(mDM) + "MeV.txt";
		// 	if(libphysica::File_Exists(file_path) == false)
		// 	{
		// 		std::vector<std::vector<double>> daily_modulation = Daily_Modulation_NREFT(25, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 		if(mpi_rank == 0)
		// 			libphysica::Export_Table(file_path, daily_modulation, {1.0, 1.0 / gram / year}, "#Time [h]\tRate [1/gr/year]");
		// 	}
		// }

		// // Tabulate the response function
		// auto l_list = libphysica::Linear_Space(0.01 * keV, 25 * keV, 100);
		// std::vector<std::vector<double>> response_function;

		// // 1. lPerpendicular
		// for(auto& l : l_list)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "l = " << l / keV << " keV" << std::endl;
		// 	std::vector<double> row = {l};
		// 	Eigen::Vector3d lVec({0, 0, l});
		// 	double W = 0.0;
		// 	for(int band = 0; band < 4; band++)
		// 	{

		// 		double w_band = graphene.Material_Response_Function(band, lVec);
		// 		W += w_band;
		// 		row.push_back(w_band);
		// 	}
		// 	row.push_back(W);
		// 	response_function.push_back(row);
		// }
		// if(mpi_rank == 0)
		// 	libphysica::Export_Table(results_path + "Response_Function_lPerp_" + cfg.carbon_wavefunctions + ".txt", response_function, {keV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV}, "#l [keV]\tW_pi [eV^-3]\tW_sigma1 [eV^-3]\tW_sigma2 [eV^-3]\tW_sigma3 [eV^-3]\tW_tot [eV^-3]");

		// // 2. lParallel (average)
		// response_function.clear();
		// for(auto& l : l_list)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "l = " << l / keV << " keV" << std::endl;
		// 	std::vector<double> row = {l};
		// 	double W				= 0.0;
		// 	for(int band = 0; band < 4; band++)
		// 	{
		// 		std::function<double(double)> integrand = [&graphene, l, band](double phi) {
		// 			Eigen::Vector3d lVec = Spherical_Coordinates(l, M_PI / 2.0, phi);
		// 			return graphene.Material_Response_Function(band, lVec);
		// 		};
		// 		double W_band = 1.0 / 2.0 / M_PI * libphysica::Integrate(integrand, 0.0, 2 * M_PI, "Gauss-Kronrod");
		// 		W += W_band;
		// 		row.push_back(W_band);
		// 	}
		// 	row.push_back(W);
		// 	response_function.push_back(row);
		// }
		// if(mpi_rank == 0)
		// 	libphysica::Export_Table(results_path + "Response_Function_lParallel_" + cfg.carbon_wavefunctions + ".txt", response_function, {keV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV}, "#l [keV]\tW_pi [eV^-3]\tW_sigma1 [eV^-3]\tW_sigma2 [eV^-3]\tW_sigma3 [eV^-3]\tW_tot [eV^-3]");

		// // 2. lNorm
		// std::ofstream f;
		// if(mpi_rank == 0)
		// 	f.open(results_path + "Response_Function_lNorm_" + cfg.carbon_wavefunctions + ".txt");
		// response_function.clear();
		// for(auto& l : l_list)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "l = " << l / keV << " keV" << std::endl;
		// 	std::vector<double> row = {l};
		// 	double W				= 0.0;
		// 	for(int band = 0; band < 4; band++)
		// 	{
		// 		std::function<double(double, double)> integrand = [&graphene, l, band](double cos_theta, double phi) {
		// 			Eigen::Vector3d lVec = Spherical_Coordinates(l, acos(cos_theta), phi);
		// 			return graphene.Material_Response_Function(band, lVec);
		// 		};
		// 		double W_band = l * l * libphysica::Integrate_2D(integrand, -1.0, 1.0, 0.0, 2 * M_PI, "Gauss-Legendre", 250);
		// 		W += W_band;
		// 		row.push_back(W_band);
		// 	}
		// 	row.push_back(W);
		// 	if(mpi_rank == 0)
		// 		f << row[0] / keV << "\t" << row[1] * eV << "\t" << row[2] * eV << "\t" << row[3] * eV << "\t" << row[4] * eV << "\t" << row[5] * eV << std::endl;
		// 	response_function.push_back(row);
		// }
		// if(mpi_rank == 0)
		// 	f.close();

		// std::vector<std::vector<double>> interpolation_list;
		// for(auto& row : response_function)
		// 	interpolation_list.push_back({row[0], row[5]});
		// libphysica::Interpolation W(interpolation_list);
		// if(mpi_rank == 0)
		// 	std::cout << "Integral = " << W.Integrate(l_list[0], l_list[l_list.size() - 1]) << std::endl;
		// if(mpi_rank == 0)
		// 	libphysica::Export_Table(results_path + "Response_Function_lNorm_" + cfg.carbon_wavefunctions + ".txt", response_function, {keV, 1.0 / eV, 1.0 / eV, 1.0 / eV, 1.0 / eV, 1.0 / eV}, "#l [keV]\tW_pi [eV^-1]\tW_sigma1 [eV^-1]\tW_sigma2 [eV^-1]\tW_sigma3 [eV^-1]\tW_tot [eV^-1]");

		// Tabulate W as a function of lx and ly.
		double lMax				   = graphene.b;
		std::vector<double> l_list = libphysica::Linear_Space(-lMax, lMax, 400);
		std::vector<std::vector<double>> response_function;
		if(mpi_rank == 0)
			std::cout << "Tabulating W(lx, ly)..." << std::endl;
		for(auto& lx : l_list)
			for(auto& ly : l_list)
			{
				Eigen::Vector3d lVec = {lx, ly, 0.0};
				double W			 = graphene.Material_Response_Function(lVec);
				response_function.push_back({lx, ly, W});
			}
		if(mpi_rank == 0)
			libphysica::Export_Table(results_path + "Response_Function_lx_ly_" + cfg.carbon_wavefunctions + ".txt", response_function, {keV, keV, 1.0 / eV / eV / eV}, "#lx [keV]\tly [eV]\tW [eV^-3]");
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}