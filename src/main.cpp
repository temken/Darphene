#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>
#include <random>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"
#include "libphysica/Special_Functions.hpp"
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
using namespace Darphene;

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
				  << LOGO << std::endl
				  << "MPI processes:\t" << mpi_processes << std::endl;
	}
	if(argc != 2)
	{
		if(mpi_rank == 0)
			std::cout << libphysica::Formatted_String("Error", "Red", true) << " in Darphene: A config file is required.\n\tCorrect usages:" << std::endl
					  << "\t>" << argv[0] << " <config_file>" << std::endl
					  << "\tor" << std::endl
					  << "\t>mpirun -n 2 " << argv[0] << " <config_file>" << std::endl
					  << std::endl;
		std::exit(EXIT_FAILURE);
	}
	////////////////////////////////////////////////////////////////////////

	// Read configuration file
	Configuration cfg(argv[1]);
	Graphene graphene(cfg.carbon_wavefunctions, cfg.graphene_work_function);
	cfg.Print_Summary(mpi_rank);
	std::string results_path = TOP_LEVEL_DIR "results/" + cfg.ID + "/";
	double rate_unit		 = 1.0 / gram / year;

	MPI_Barrier(MPI_COMM_WORLD);

	// Run modes:
	if(cfg.run_modus == "Energy-Spectrum")
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
	else if(cfg.run_modus == "Total-Rate")
	{
		if(mpi_rank == 0)
			std::cout << "\nCompute total rate R:" << std::endl;
		double R = R_Total_NREFT(*cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		std::cout << "R = " << In_Units(R, rate_unit) << " / gr / year" << std::endl;
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Total_Rate.txt";
			std::ofstream f(file_path);
			f << "R = " << In_Units(R, rate_unit) << " /gr / year" << std::endl;
			f.close();
			std::cout << "\nDone. Total rate saved to " << file_path << "." << std::endl;
		}
	}
	else if(cfg.run_modus == "Directional-Spectrum")
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
	else if(cfg.run_modus == "Daily-Modulation")
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
	else if(cfg.run_modus == "Graphene")
	{
		// 1. Graphene band structure
		if(mpi_rank == 0)
			std::cout << "\n1. Tabulate graphene band structure:" << std::endl;
		auto band_structure = graphene.Energy_Bands(cfg.grid_points);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Graphene_Band_Structure.txt";
			libphysica::Export_Table(file_path, band_structure, {keV, eV, eV, eV, eV, eV, eV, eV, eV}, "# k [keV]\tE_pi1 [eV]\tE_pi2 [eV]\tE_sigma11 [eV]\tE_sigma12 [eV]\tE_sigma21 [eV]\tE_sigma22 [eV]\tE_sigma31 [eV]\tE_sigma32 [eV]");
			std::cout << "\nDone. Tabulated band structure saved to \n\t" << file_path << "." << std::endl;
		}

		// 2. Tabulate the graphene response function
		if(mpi_rank == 0)
			std::cout << "\n2. Tabulate graphene response function:" << std::endl;
		auto l_list = libphysica::Linear_Space(0.01 * keV, 10 * keV, cfg.grid_points);
		std::vector<std::vector<double>> response_function;

		// 2.1. lPerpendicular
		if(mpi_rank == 0)
			std::cout << "\n2.1. lPerpendicular:" << std::endl;
		for(auto& l : l_list)
		{
			std::vector<double> row = {l};
			Eigen::Vector3d lVec({0, 0, l});
			double W = 0.0;
			for(int band = 0; band < 4; band++)
			{

				double w_band = graphene.Material_Response_Function(band, lVec);
				W += w_band;
				row.push_back(w_band);
			}
			row.push_back(W);
			response_function.push_back(row);
			libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size(), mpi_rank, 86, 0.0, "Blue");
		}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_lPerp.txt";
			libphysica::Export_Table(file_path, response_function, {keV, 1.0 / keV / keV / keV, 1.0 / keV / keV / keV, 1.0 / keV / keV / keV, 1.0 / keV / keV / keV, 1.0 / keV / keV / keV}, "#l [keV]\tW_pi [keV^-3]\tW_sigma1 [keV^-3]\tW_sigma2 [keV^-3]\tW_sigma3 [keV^-3]\tW_tot [keV^-3]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}

		// 2.2 lParallel (average)
		if(mpi_rank == 0)
			std::cout << "\n2.2. lParallel:" << std::endl;
		response_function.clear();
		for(auto& l : l_list)
		{
			std::vector<double> row = {l};
			double W				= 0.0;
			for(int band = 0; band < 4; band++)
			{
				std::function<double(double)> integrand = [&graphene, l, band](double phi) {
					Eigen::Vector3d lVec = Spherical_Coordinates(l, M_PI / 2.0, phi);
					return graphene.Material_Response_Function(band, lVec);
				};
				double W_band = 1.0 / 2.0 / M_PI * libphysica::Integrate(integrand, 0.0, 2 * M_PI, "Gauss-Kronrod");
				W += W_band;
				row.push_back(W_band);
			}
			row.push_back(W);
			response_function.push_back(row);
			libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size(), mpi_rank, 86, 0.0, "Blue");
		}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_lParallel.txt";
			libphysica::Export_Table(file_path, response_function, {keV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV}, "#l [keV]\tW_pi [eV^-3]\tW_sigma1 [eV^-3]\tW_sigma2 [eV^-3]\tW_sigma3 [eV^-3]\tW_tot [eV^-3]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}

		// 2.3 lNorm
		if(mpi_rank == 0)
			std::cout << "\n2.3. lNorm:" << std::endl;
		response_function.clear();
		for(auto& l : l_list)
		{
			std::vector<double> row = {l};
			double W				= 0.0;
			for(int band = 0; band < 4; band++)
			{
				std::function<double(double, double)> integrand = [&graphene, l, band](double cos_theta, double phi) {
					Eigen::Vector3d lVec = Spherical_Coordinates(l, acos(cos_theta), phi);
					return graphene.Material_Response_Function(band, lVec);
				};
				double W_band = l * l * libphysica::Integrate_2D(integrand, -1.0, 1.0, 0.0, 2 * M_PI, "Vegas", 2000);
				W += W_band;
				row.push_back(W_band);
			}
			row.push_back(W);
			response_function.push_back(row);
			libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size(), mpi_rank, 86, 0.0, "Blue");
		}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_lNorm.txt";
			libphysica::Export_Table(file_path, response_function, {keV, 1.0 / keV, 1.0 / keV, 1.0 / keV, 1.0 / keV, 1.0 / keV}, "#l [keV]\tW_pi [keV^-1]\tW_sigma1 [keV^-1]\tW_sigma2 [keV^-1]\tW_sigma1 [keV^-1]\tW_tot [keV^-1]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}

		// 2.4 Tabulate W as a function of lx,ly
		if(mpi_rank == 0)
			std::cout << "\n2.4. W(lx,ly,lz=0):" << std::endl;
		double lMax = graphene.b;
		double lz	= 0.0 * eV;
		l_list		= libphysica::Linear_Space(-lMax, lMax, cfg.grid_points);
		response_function.clear();
		for(auto& lx : l_list)
			for(auto& ly : l_list)
			{
				Eigen::Vector3d lVec = {lx, ly, lz};
				double W_pi			 = graphene.Material_Response_Function(0, lVec);
				double W_s1			 = graphene.Material_Response_Function(1, lVec);
				double W_s2			 = graphene.Material_Response_Function(2, lVec);
				double W_s3			 = graphene.Material_Response_Function(3, lVec);
				double W_tot		 = W_pi + W_s1 + W_s2 + W_s3;
				response_function.push_back({lx, ly, W_pi, W_s1, W_s2, W_s3, W_tot});
				libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size() / l_list.size(), mpi_rank, 86, 0.0, "Blue");
			}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_lx_ly_lz=0eV.txt";
			libphysica::Export_Table(file_path, response_function, {keV, keV, std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3)}, "#lx [keV]\tly [keV]\tW_pi [keV^-3]\tW_s1 [keV^-3]\tW_s2 [keV^-3]\tW_s3 [keV^-3]\tW [keV^-3]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}

		// 2.5 Tabulate W as a function of lx,ly, lz.
		if(mpi_rank == 0)
			std::cout << "\n2.5. W(lx,ly,lz):" << std::endl;
		response_function.clear();
		int counter = 0;
		for(auto& lx : l_list)
			for(auto& ly : l_list)
				for(auto& lz : l_list)
				{
					Eigen::Vector3d lVec = {lx, ly, lz};
					double W_pi			 = graphene.Material_Response_Function(0, lVec);
					double W_s1			 = graphene.Material_Response_Function(1, lVec);
					double W_s2			 = graphene.Material_Response_Function(2, lVec);
					double W_s3			 = graphene.Material_Response_Function(3, lVec);
					double W_tot		 = W_pi + W_s1 + W_s2 + W_s3;
					response_function.push_back({lx, ly, lz, W_pi, W_s1, W_s2, W_s3, W_tot});
					if(counter++ % 10 == 0)
						libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size() / l_list.size() / l_list.size(), mpi_rank, 86, 0.0, "Blue");
				}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_lx_ly_lz.txt";
			libphysica::Export_Table(file_path, response_function, {keV, keV, keV, std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3)}, "#lx [keV]\tly [keV]\tlz [keV]\tW_pi [keV^-3]\tW_s1 [keV^-3]\tW_s2 [keV^-3]\tW_s3 [keV^-3]\tW [keV^-3]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}

		// 2.6 Tabulate W as a function of l and theta and integrate over phi.
		if(mpi_rank == 0)
			std::cout << "\n2.6. W(l,theta):" << std::endl;
		l_list			= libphysica::Linear_Space(0.01 * keV, 10 * keV, cfg.grid_points);
		auto theta_list = libphysica::Linear_Space(0.0, 0.5 * M_PI, cfg.grid_points);
		response_function.clear();
		for(auto& l : l_list)
		{
			for(auto& theta : theta_list)
			{
				std::function<double(double)> integrand = [&graphene, l, theta](double phi) {
					Eigen::Vector3d lVec = Spherical_Coordinates(l, theta, phi);
					return graphene.Material_Response_Function(lVec);
				};
				double W = l * l * libphysica::Integrate(integrand, 0.0, 2 * M_PI, "Gauss-Kronrod");
				response_function.push_back({l, theta, W});
			}
			libphysica::Print_Progress_Bar(1.0 * response_function.size() / l_list.size() / theta_list.size(), mpi_rank, 86, 0.0, "Blue");
		}
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Response_Function_l_theta.txt";
			libphysica::Export_Table(file_path, response_function, {keV, deg, std::pow(keV, -1)}, "#l [keV]\ttheta [deg]\tW [keV^-1]");
			std::cout << "\nDone. Tabulated response function saved to \n\t" << file_path << "." << std::endl;
		}
	}
	else if(cfg.run_modus == "Exclusion-Limit")
	{
		if(mpi_rank == 0)
			std::cout << "\nCompute exclusion limit for standard SI interactions:" << std::endl;
		DM_Detector_Graphene detector("Graphene-Detector", cfg.exposure, graphene);
		auto DM_masses		 = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);
		auto exclusion_limit = detector.Upper_Limit_Curve(*cfg.DM, *cfg.DM_distr, DM_masses, cfg.constraints_certainty);
		if(mpi_rank == 0)
		{
			std::string file_path = results_path + "Exclusion_Limit.txt";
			libphysica::Export_Table(file_path, exclusion_limit, {MeV, cm * cm}, "#m_DM [MeV]\t#sigma_e_SI [cm^2]");
			std::cout << "\nDone. Exclusion limit saved to \n\t" << file_path << "." << std::endl;
		}
	}
	else if(cfg.run_modus == "Custom")
	{
	}
	else
	{
		throw std::runtime_error("Error: Run modus " + cfg.run_modus + " not recognized.");
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << libphysica::Formatted_String("\n[Finished in " + libphysica::Time_Display(durationTotal) + "]\a", "Green", true) << std::endl;
	MPI_Finalize();
	return 0;
}
