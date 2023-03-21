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
using namespace DarPhene;

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
			std::cout << libphysica::Formatted_String("Error", "Red", true) << " in DarPhene: A config file is required.\n\tCorrect usages:" << std::endl
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
		// if(mpi_rank == 0)
		// 	std::cout << mpi_rank << "\t" << In_Units(R_Total_NREFT(*cfg.DM_NREFT, *cfg.DM_distr, graphene, 0, cfg.MC_points), rate_unit) << std::endl;
		// else
		// 	std::cout << mpi_rank << "\t" << In_Units(R_Total_NREFT_alt(*cfg.DM_NREFT, *cfg.DM_distr, graphene, 0, cfg.MC_points), rate_unit) << std::endl;

		// for(auto overlap : {"s", "sPrime", "Sss", "Ssigma", "Ssp"})
		// 	std::cout << overlap << "\t" << graphene.Overlap_Integral(overlap) << std::endl;

		// // 1. Band structure
		// auto band_structure	  = graphene.Energy_Bands(200);
		// std::string file_path = results_path + "Graphene_Band_Structure.txt";
		// libphysica::Export_Table(file_path, band_structure, {keV, eV, eV, eV, eV, eV, eV, eV, eV}, "# k [keV]\tE_pi1 [eV]\tE_pi2 [eV]\tE_sigma11 [eV]\tE_sigma12 [eV]\tE_sigma21 [eV]\tE_sigma22 [eV]\tE_sigma31 [eV]\tE_sigma32 [eV]");

		// 2. Tabulate the response function
		// auto l_list = libphysica::Linear_Space(0.01 * keV, 10 * keV, 100);
		// std::vector<std::vector<double>> response_function;

		// // 2.1. lPerpendicular
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
		// 	libphysica::Export_Table(results_path + "Response_Function_lPerp.txt", response_function, {keV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV}, "#l [keV]\tW_pi [eV^-3]\tW_sigma1 [eV^-3]\tW_sigma2 [eV^-3]\tW_sigma3 [eV^-3]\tW_tot [eV^-3]");

		// // 2.2 lParallel (average)
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
		// 	libphysica::Export_Table(results_path + "Response_Function_lParallel.txt", response_function, {keV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV, 1.0 / eV / eV / eV}, "#l [keV]\tW_pi [eV^-3]\tW_sigma1 [eV^-3]\tW_sigma2 [eV^-3]\tW_sigma3 [eV^-3]\tW_tot [eV^-3]");

		// // 2.3 lNorm
		// std::ofstream f;
		// if(mpi_rank == 0)
		// 	f.open(results_path + "Response_Function_lNorm_total_Vegas.txt");
		// response_function.clear();
		// for(auto& l : l_list)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "l = " << l / keV << " keV" << std::endl;
		// 	std::vector<double> row = {l};
		// 	double W				= 0.0;
		// 	// for(int band = 0; band < 4; band++)
		// 	// {
		// 	std::function<double(double, double)> integrand = [&graphene, l](double cos_theta, double phi) {
		// 		Eigen::Vector3d lVec = Spherical_Coordinates(l, acos(cos_theta), phi);
		// 		// return graphene.Material_Response_Function(band, lVec);
		// 		return graphene.Material_Response_Function(lVec);
		// 	};
		// 	double W_band = l * l * libphysica::Integrate_2D(integrand, -1.0, 1.0, 0.0, 2 * M_PI, "Vegas", 2000);
		// 	W += W_band;
		// 	// row.push_back(W_band);
		// 	// }
		// 	row.push_back(W);
		// 	if(mpi_rank == 0)
		// 		f << row[0] / keV << "\t" << row[1] * eV << std::endl;	 //<< "\t" << row[2] * eV << "\t" << row[3] * eV << "\t" << row[4] * eV << "\t" << row[5] * eV << std::endl;
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
		// 	libphysica::Export_Table(results_path + "Response_Function_lNorm.txt", response_function, {keV, 1.0 / eV, 1.0 / eV, 1.0 / eV, 1.0 / eV, 1.0 / eV}, "#l [keV]\tW_pi [eV^-1]\tW_sigma1 [eV^-1]\tW_sigma2 [eV^-1]\tW_sigma3 [eV^-1]\tW_tot [eV^-1]");

		// // 2.4 Tabulate W as a function of lx,ly
		// double lMax = graphene.b;
		// // std::vector<double> l_list = libphysica::Linear_Space(-lMax, lMax, 200);
		// // std::vector<std::vector<double>> response_function;
		// l_list = libphysica::Linear_Space(-lMax, lMax, 200);
		// response_function.clear();
		// if(mpi_rank == 0)
		// 	std::cout << "Tabulating W(lx, ly)..." << std::endl;
		// for(auto& lx : l_list)
		// 	for(auto& ly : l_list)
		// 	{
		// 		Eigen::Vector3d lVec = {lx, ly, 91.0 * eV};
		// 		double W_pi			 = graphene.Material_Response_Function(0, lVec);
		// 		double W_s1			 = graphene.Material_Response_Function(1, lVec);
		// 		double W_s2			 = graphene.Material_Response_Function(2, lVec);
		// 		double W_s3			 = graphene.Material_Response_Function(3, lVec);
		// 		double W_tot		 = W_pi + W_s1 + W_s2 + W_s3;
		// 		response_function.push_back({lx, ly, W_pi, W_s1, W_s2, W_s3, W_tot});
		// 	}
		// if(mpi_rank == 0)
		// 	libphysica::Export_Table(results_path + "Response_Function_lx_ly_lz=91eV.txt", response_function, {keV, keV, std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3)}, "#lx [keV]\tly [keV]\tW_pi [keV^-3]\tW_s1 [keV^-3]\tW_s2 [keV^-3]\tW_s3 [keV^-3]\tW [keV^-3]");

		// // 2.5 Tabulate W as a function of lx,ly, lz.
		// double lMax = graphene.b;
		// // std::vector<double>
		// l_list = libphysica::Linear_Space(-lMax, lMax, 50);
		// response_function.clear();
		// if(mpi_rank == 0)
		// 	std::cout << "Tabulating W(lx, ly, lz)..." << std::endl;
		// for(auto& lx : l_list)
		// 	for(auto& ly : l_list)
		// 		for(auto& lz : l_list)
		// 		{
		// 			Eigen::Vector3d lVec = {lx, ly, lz};
		// 			double W_pi			 = graphene.Material_Response_Function(0, lVec);
		// 			double W_s1			 = graphene.Material_Response_Function(1, lVec);
		// 			double W_s2			 = graphene.Material_Response_Function(2, lVec);
		// 			double W_s3			 = graphene.Material_Response_Function(3, lVec);
		// 			double W_tot		 = W_pi + W_s1 + W_s2 + W_s3;
		// 			response_function.push_back({lx, ly, lz, W_pi, W_s1, W_s2, W_s3, W_tot});
		// 		}
		// if(mpi_rank == 0)
		// 	libphysica::Export_Table(results_path + "Response_Function_lx_ly_lz.txt", response_function, {keV, keV, keV, std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3), std::pow(keV, -3)}, "#lx [keV]\tly [keV]\tlz [keV]\tW_pi [keV^-3]\tW_s1 [keV^-3]\tW_s2 [keV^-3]\tW_s3 [keV^-3]\tW [keV^-3]");

		// // 2.6 Tabulate W as a function of l and theta and integrate over phi.
		// l_list			= libphysica::Linear_Space(1.0 * keV, 10 * keV, 150);
		// auto theta_list = libphysica::Linear_Space(0.0, 0.5 * M_PI, 150);
		// std::ofstream f;
		// if(mpi_rank == 0)
		// 	f.open(results_path + "Response_Function_l_theta_5eV.txt");
		// response_function.clear();
		// for(auto& l : l_list)
		// {
		// 	if(mpi_rank == 0)
		// 		std::cout << "l = " << l / keV << " keV" << std::endl;
		// 	for(auto& theta : theta_list)
		// 	{

		// 		std::vector<double> row = {l, theta};
		// 		double W				= 0.0;
		// 		// for(int band = 0; band < 4; band++)
		// 		// {
		// 		std::function<double(double)> integrand = [&graphene, l, theta](double phi) {
		// 			Eigen::Vector3d lVec = Spherical_Coordinates(l, theta, phi);
		// 			// return graphene.Material_Response_Function(band, lVec);
		// 			return graphene.Material_Response_Function(lVec);
		// 		};
		// 		double W_band = l * l * libphysica::Integrate(integrand, 0.0, 2 * M_PI, "Gauss-Kronrod");
		// 		W += W_band;
		// 		// row.push_back(W_band);
		// 		// }
		// 		row.push_back(W);
		// 		if(mpi_rank == 0)
		// 			// f << row[0] / keV << "\t" << row[1] / deg << "\t" << row[2] * eV << "\t" << row[3] * eV << "\t" << row[4] * eV << "\t" << row[5] * eV << "\t" << row[6] * eV << std::endl;
		// 			f << row[0] / keV << "\t" << row[1] / deg << "\t" << row[2] * eV << std::endl;
		// 		response_function.push_back(row);
		// 	}
		// }
		// if(mpi_rank == 0)
		// 	f.close();

		// // 1. DM model and mass
		// // std::vector<double> DM_masses = libphysica::Log_Space(mMin, 100.0 * MeV, 10);
		// std::vector<int> DM_masses = {100};	  // {2, 5, 10, 20, 50, 100};

		// std::vector<int> operators			  = {1};	 //= {1, 3};
		// std::vector<std::string> interactions = {"L"};	 //{"C", "L"};
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
		// 			std::string file_path = results_path + "Daily_Modulation_O" + std::to_string(op) + inter + "_m_" + std::to_string(mDM) + "_MeV_new.txt";
		// 			// if(libphysica::File_Exists(file_path) == false)
		// 			// {
		// 			std::vector<std::vector<double>> daily_modulation = Daily_Modulation_NREFT(cfg.grid_points, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 			if(mpi_rank == 0)
		// 				libphysica::Export_Table(file_path, daily_modulation, {1.0, 1.0 / gram / year}, "#Time [h]\tRate [1/gr/year]");
		// 			// }
		// 		}

		// // std::vector<double> DM_masses		  = libphysica::Log_Space(2.0 * MeV, 100.0 * MeV, 10);
		// std::vector<double> DM_masses		  = {1.425229299235770242e+00 * MeV, 1.658889390466721236e+00 * MeV, 2.000000000000000000e+00 * MeV, 2.247412567913343651e+00 * MeV, 2.615865999184929347e+00 * MeV, 3.044725754134705475e+00 * MeV, 3.543895184531499432e+00 * MeV, 4.124901253221343822e+00 * MeV, 5.000000000000000000e+00 * MeV, 5.588289955415879362e+00 * MeV, 6.504465578178933782e+00 * MeV, 7.570844175097216855e+00 * MeV, 8.812050864853802068e+00 * MeV, 1.000000000000000000e+01 * MeV, 1.193829677649106813e+01 * MeV, 1.389552816054473894e+01 * MeV, 1.617363904378025197e+01 * MeV, 2.000000000000000000e+01 * MeV, 2.191155110447665066e+01 * MeV, 2.550385378365185574e+01 * MeV, 2.968509872790354365e+01 * MeV, 3.455184043794348980e+01 * MeV, 4.021646310129752067e+01 * MeV, 5.000000000000000000e+01 * MeV, 5.448403319129700151e+01 * MeV, 6.341645141546277387e+01 * MeV, 7.381330042160217886e+01 * MeV, 8.591466721206701607e+01 * MeV, 1.000000000000000000e+02 * MeV};
		// std::vector<int> operators			  = {3};	 //= {1, 3};
		// std::vector<std::string> interactions = {"C"};	 //{"C", "L"};
		// // int tot								  = DM_masses.size() * operators.size() * interactions.size();
		// // int i								  = 0;
		// for(auto& op : operators)
		// 	for(auto& inter : interactions)
		// 	{
		// 		std::ofstream f;
		// 		f.open(results_path + "Average_Rate_Mass_O" + std::to_string(op) + inter + ".txt", std::ios_base::app);
		// 		if(mpi_rank == 0)
		// 			f << std::endl;
		// 		// for(int j = DM_masses.size() - 1; j >= 0; j--)
		// 		for(int j = 0; j >= 0; j--)
		// 		{
		// 			double mDM = DM_masses[j];
		// 			if(mpi_rank == 0)
		// 				std::cout << j << ")\tmDM = " << mDM / MeV << " MeV" << std::endl;
		// 			DM_Particle_NREFT DM(mDM);
		// 			DM.Set_Coupling(op, 1.0, inter);
		// 			std::vector<std::vector<double>> daily_nreft = Daily_Modulation_NREFT(cfg.grid_points, DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 			double average								 = 0.0;
		// 			for(int k = 0; k < daily_nreft.size(); k++)
		// 				average += daily_nreft[k][1];
		// 			average /= daily_nreft.size();
		// 			if(mpi_rank == 0)
		// 			{
		// 				std::cout << "\nAverage rate = " << average / rate_unit << " /gr/year" << std::endl;
		// 				f << mDM / MeV << "\t" << average / rate_unit << std::endl;
		// 				std::string file_path = results_path + "Daily_Modulation_O" + std::to_string(op) + inter + "_m=" + std::to_string(libphysica::Round(mDM / MeV)) + ".txt";
		// 				libphysica::Export_Table(file_path, daily_nreft, {1.0, rate_unit}, "# t[h]\tR[1/gr/year]");
		// 				std::cout << "\nDone. Table saved to " << file_path << "." << std::endl;
		// 			}
		// 			// double rate = R_Average_NREFT(DM, *cfg.DM_distr, graphene, cfg.MC_points);
		// 			// f << mDM / MeV << "\t" << In_Units(rate, rate_unit) << std::endl;
		// 			// std::cout << ++i << "/" << tot << "\t" << op << "\t" << inter << "\t" << mDM / MeV << "\t" << In_Units(rate, rate_unit) << std::endl;
		// 		}
		// 		f.close();
		// 	}
	}
	else
	{
		// throw std::runtime_error("Error: Run modus " + cfg.run_modus + " not recognized.");
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