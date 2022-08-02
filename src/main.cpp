#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <random>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

#include "graphene/Configuration.hpp"
#include "graphene/DM_Particle_NREFT.hpp"
#include "graphene/Direct_Detection_Standard.hpp"
#include "graphene/Graphene.hpp"
#include "graphene/Hydrogenic_Wavefunctions.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace graphene;

int main(int argc, char* argv[])
{
	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << std::endl;
	////////////////////////////////////////////////////////////////////////

	Configuration cfg(argv[1]);
	Graphene graphene(cfg.graphene_work_function);
	cfg.Print_Summary();
	std::string results_path = TOP_LEVEL_DIR "results/" + cfg.ID + "/";
	double rate_unit		 = 1.0 / gram / year;

	if(cfg.run_modus == "Energy-Spectrum" || cfg.run_modus == "All")
	{
		std::cout << "\nTabulate dR/dlnE:" << std::endl;
		auto spectrum_nreft	  = Tabulate_dR_dlnE_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "dR_dlnE.txt";
		libphysica::Export_Table(file_path, spectrum_nreft, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});
		std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
	}
	if(cfg.run_modus == "Total-Rate" || cfg.run_modus == "All")
	{
		std::cout << "\nCompute total rate R:" << std::endl;
		double R			  = R_Total_NREFT(*cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		std::string file_path = results_path + "Total_Rate.txt";
		std::ofstream f(file_path);
		f << "R = " << In_Units(R, rate_unit) << " /gr / year" << std::endl;
		std::cout << "R = " << In_Units(R, rate_unit) << " / gr / year" << std::endl;
		f.close();
		std::cout << "\nDone. Total rate saved to " << file_path << "." << std::endl;
	}
	if(cfg.run_modus == "Directional-Spectrum" || cfg.run_modus == "All")
	{
		std::cout << "\nTabulate dR/(dcos dphi):" << std::endl;
		auto spectrum_nreft	  = Tabulate_dR_dcos_dphi_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "dR_dcos_dphi.txt";
		libphysica::Export_Table(file_path, spectrum_nreft, {1.0, 1.0, rate_unit});
		std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
	}
	if(cfg.run_modus == "Daily-Modulation" || cfg.run_modus == "All")
	{
		std::cout << "\nCalculate daily modulation:" << std::endl;
		auto daily_nreft	  = Daily_Modulation_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "Daily_Modulation.txt";
		libphysica::Export_Table(file_path, daily_nreft, {1.0, rate_unit});
		std::cout << "\nDone. Table saved to " << file_path << "." << std::endl;
	}
	else if(cfg.run_modus == "Custom")
	{

		auto cos_list = libphysica::Linear_Space(-1.0, 1.0, 25);
		std::ofstream f(cfg.results_path + "dR_dcos_t=0.txt");
		for(auto cos : cos_list)
		{
			double dR_dcos = dR_dcos_NREFT(cos, *cfg.DM_NREFT, *cfg.DM_distr, graphene, 0, cfg.MC_points);
			f << cos << "\t" << dR_dcos / rate_unit << std::endl;
			std::cout << cos << "\t" << dR_dcos / rate_unit << std::endl;
		}
		f.close();

		///////////////////////////////////////////////
		// double mDM				= 20.0 * MeV;
		// std::string interaction = "Long-Range";
		// std::ofstream f(TOP_LEVEL_DIR "results/Total_Rate_" + interaction + ".txt");
		// for(int op = 1; op < 16; op++)
		// {
		// 	if(op == 2)
		// 		continue;
		// 	std::cout << "\nCalculate " << interaction << " interaction for op = " << op << std::endl;
		// 	DM_Particle_NREFT DM_nreft(mDM);
		// 	DM_nreft.Set_Coupling(op, 1.0, interaction);

		// 	std::cout << "\nTabulate dR/dlnE:" << std::endl;
		// 	auto spectrum_nreft	  = Tabulate_dR_dlnE_NREFT(cfg.grid_points, DM_nreft, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		// 	std::string file_path = TOP_LEVEL_DIR "results/all_spectra/O" + std::to_string(op) + "_dR_dlnE_" + interaction + ".txt";
		// 	libphysica::Export_Table(file_path, spectrum_nreft, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});
		// 	std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;

		// 	std::cout << "\nCompute total rate R:" << std::endl;
		// 	double R = R_Total_NREFT(DM_nreft, *cfg.DM_distr, graphene, cfg.MC_points);
		// 	f << "O" << op << ":\tR = " << In_Units(R, rate_unit) << " /gr / year" << std::endl;
		// 	std::cout << "O" << op << ":\tR = " << In_Units(R, rate_unit) << " /gr / year" << std::endl;
		// 	std::cout << "\nDone. Total rate saved to " << file_path << "." << std::endl;
		// }
		// f.close();
		///////////////////////////////////////////////

		// extern double R_Total_NREFT(DM_Particle_NREFT& DM, obscura::DM_Distribution& DM_distr, Graphene& graphene, unsigned int MC_points);

		// double R = R_Total_NREFT(*cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points);
		// std::cout << "R = " << In_Units(R, 1.0 / gram / year) << std::endl;

		// obscura::DM_Particle_SI DM_standard(100 * MeV);
		// DM_standard.Set_Sigma_Electron(1.0e-36 * cm * cm);
		// DM_Particle_NREFT DM_nreft = DM_Dark_Photon(100.0 * MeV, pb, "Contact");

		// int band = 0;
		// std::cout << "\nQty\tStd\tSimple\tNREFT" << std::endl;
		// double R_standard		 = R_Total_Standard(DM_standard, *cfg.DM_distr, graphene, band, 1e5);		   //, band);
		// double R_standard_simple = R_Total_Standard_Simple(DM_standard, *cfg.DM_distr, graphene, band, 1e4);   //, band);
		// double R_nreft			 = R_Total_NREFT(DM_nreft, *cfg.DM_distr, graphene, 0, 1e7);				   //, band);

		// std::cout << "R\t" << In_Units(R_standard, rate_unit) << "\t" << In_Units(R_standard_simple, rate_unit) << "\t" << In_Units(R_nreft, rate_unit) << std::endl;

		// double Ee					   = 5.0 * eV;
		// double dR_dlnE_standard		   = dR_dlnE_Standard(Ee, DM_standard, SHM, graphene, band, 1e4);		   //, band);
		// double dR_dlnE_standard_simple = dR_dlnE_Standard_Simple(Ee, DM_standard, SHM, graphene, band, 1e4);   //, band);
		// double dR_dlnE_nreft		   = dR_dlnE_NREFT(Ee, DM_nreft, SHM, graphene, 0, 1e5);				   //, band);

		// std::cout << "dRdE\t" << In_Units(dR_dlnE_standard, rate_unit) << "\t" << In_Units(dR_dlnE_standard_simple, rate_unit) << "\t" << In_Units(dR_dlnE_nreft, rate_unit) << std::endl;
		// double cos					 = 0.3;
		// double phi					 = 0.3;
		// double dR_dcos_dphi_standard = dR_dcos_dphi_Standard(cos, phi, DM_standard, SHM, graphene, band, 1e5);	 //, band);
		// double dR_dcos_dphi_nreft	 = dR_dcos_dphi_NREFT(cos, phi, DM_nreft, SHM, graphene, 0, 1e5);			 //, band);

		// std::cout << "dRdcdp\t" << In_Units(dR_dcos_dphi_standard, rate_unit) << "\t\t" << In_Units(dR_dcos_dphi_nreft, rate_unit) << std::endl;
		// 1. Energy spectra
		// int points = 15;
		// std::cout << "Tabulate dR/dlnE for standard" << std::endl;
		// auto spectrum_std = Tabulate_dR_dlnE_Standard(points, DM_standard, SHM, graphene, "Full", 1e4);
		// libphysica::Export_Table("dR_dlnE_Standard.txt", spectrum_std, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});
		// std::cout << "Tabulate dR/dlnE for simple" << std::endl;
		// auto spectrum_sim = Tabulate_dR_dlnE_Standard(points, DM_standard, SHM, graphene, "Simplified", 1e3);
		// libphysica::Export_Table("dR_dlnE_Simple.txt", spectrum_sim, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});
		// std::cout << "Tabulate dR/dlnE for NREFT" << std::endl;
		// auto spectrum_nreft = Tabulate_dR_dlnE_NREFT(points, DM_nreft, SHM, graphene, 1e4);
		// libphysica::Export_Table("dR_dlnE_NREFT.txt", spectrum_nreft, {eV, rate_unit, rate_unit, rate_unit, rate_unit, rate_unit});

		// // 1. Daily modulation
		// int points = 4;
		// // std::cout << "Tabulate daily modulation for standard" << std::endl;
		// // auto daily_std = Daily_Modulation_Standard(points, DM_standard, SHM, graphene, 1e5);
		// // libphysica::Export_Table("Daily_Modulation_Standard.txt", daily_std, {1.0, rate_unit});
		// std::cout << "Tabulate daily modulation for nreft" << std::endl;
		// auto daily_nreft = Daily_Modulation_NREFT(points, DM_nreft, SHM, graphene, 2e6);
		// libphysica::Export_Table("Daily_Modulation_NREFT.txt", daily_nreft, {1.0, rate_unit});
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}