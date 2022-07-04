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
#include "graphene/Direct_Detection_Graphene.hpp"
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
	double dim				 = 1.0 / kg / year;

	if(cfg.run_modus == "Energy-Spectrum" || cfg.run_modus == "All")
	{
		std::cout << "Tabulate dR/dlnE:" << std::endl;
		auto spectrum_nreft	  = Tabulate_dR_dlnE_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "dR_dlnE.txt";
		libphysica::Export_Table(file_path, spectrum_nreft, {eV, dim, dim, dim, dim, dim});
		std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
	}
	if(cfg.run_modus == "Directional-Spectrum" || cfg.run_modus == "All")
	{
		std::cout << "Tabulate dR/(dcos dphi):" << std::endl;
		auto spectrum_nreft	  = Tabulate_dR_dcos_dphi_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "dR_dcos_dphi.txt";
		libphysica::Export_Table(file_path, spectrum_nreft, {1.0, 1.0, dim});
		std::cout << "\nDone. Tabulated spectrum saved to " << file_path << "." << std::endl;
	}
	if(cfg.run_modus == "Daily-Modulation" || cfg.run_modus == "All")
	{
		std::cout << "Calculate daily modulation:" << std::endl;
		auto daily_nreft	  = Daily_Modulation_NREFT(cfg.grid_points, *cfg.DM_NREFT, *cfg.DM_distr, graphene, cfg.MC_points, cfg.threads);
		std::string file_path = results_path + "Daily_Modulation.txt";
		libphysica::Export_Table(file_path, daily_nreft, {1.0, dim});
		std::cout << "\nDone. Table saved to " << file_path << "." << std::endl;
	}
	else if(cfg.run_modus == "Custom")
	{
		// obscura::DM_Particle_SI DM_standard(100 * MeV);
		// DM_standard.Set_Sigma_Electron(1.0e-36 * cm * cm);
		// DM_Particle_NREFT DM_nreft = DM_Dark_Photon(100.0 * MeV, pb, "Contact");

		// int band = 0;
		// std::cout << "\nQty\tStd\tSimple\tNREFT" << std::endl;
		// double R_standard		 = R_Total_Standard(DM_standard, *cfg.DM_distr, graphene, band, 1e5);		   //, band);
		// double R_standard_simple = R_Total_Standard_Simple(DM_standard, *cfg.DM_distr, graphene, band, 1e4);   //, band);
		// double R_nreft			 = R_Total_NREFT(DM_nreft, *cfg.DM_distr, graphene, 0, 1e7);				   //, band);

		// std::cout << "R\t" << In_Units(R_standard, dim) << "\t" << In_Units(R_standard_simple, dim) << "\t" << In_Units(R_nreft, dim) << std::endl;

		// double Ee					   = 5.0 * eV;
		// double dR_dlnE_standard		   = dR_dlnE_Standard(Ee, DM_standard, SHM, graphene, band, 1e4);		   //, band);
		// double dR_dlnE_standard_simple = dR_dlnE_Standard_Simple(Ee, DM_standard, SHM, graphene, band, 1e4);   //, band);
		// double dR_dlnE_nreft		   = dR_dlnE_NREFT(Ee, DM_nreft, SHM, graphene, 0, 1e5);				   //, band);

		// std::cout << "dRdE\t" << In_Units(dR_dlnE_standard, dim) << "\t" << In_Units(dR_dlnE_standard_simple, dim) << "\t" << In_Units(dR_dlnE_nreft, dim) << std::endl;
		// double cos					 = 0.3;
		// double phi					 = 0.3;
		// double dR_dcos_dphi_standard = dR_dcos_dphi_Standard(cos, phi, DM_standard, SHM, graphene, band, 1e5);	 //, band);
		// double dR_dcos_dphi_nreft	 = dR_dcos_dphi_NREFT(cos, phi, DM_nreft, SHM, graphene, 0, 1e5);			 //, band);

		// std::cout << "dRdcdp\t" << In_Units(dR_dcos_dphi_standard, dim) << "\t\t" << In_Units(dR_dcos_dphi_nreft, dim) << std::endl;
		// 1. Energy spectra
		// int points = 15;
		// std::cout << "Tabulate dR/dlnE for standard" << std::endl;
		// auto spectrum_std = Tabulate_dR_dlnE_Standard(points, DM_standard, SHM, graphene, "Full", 1e4);
		// libphysica::Export_Table("dR_dlnE_Standard.txt", spectrum_std, {eV, dim, dim, dim, dim, dim});
		// std::cout << "Tabulate dR/dlnE for simple" << std::endl;
		// auto spectrum_sim = Tabulate_dR_dlnE_Standard(points, DM_standard, SHM, graphene, "Simplified", 1e3);
		// libphysica::Export_Table("dR_dlnE_Simple.txt", spectrum_sim, {eV, dim, dim, dim, dim, dim});
		// std::cout << "Tabulate dR/dlnE for NREFT" << std::endl;
		// auto spectrum_nreft = Tabulate_dR_dlnE_NREFT(points, DM_nreft, SHM, graphene, 1e4);
		// libphysica::Export_Table("dR_dlnE_NREFT.txt", spectrum_nreft, {eV, dim, dim, dim, dim, dim});

		// // 1. Daily modulation
		// int points = 4;
		// // std::cout << "Tabulate daily modulation for standard" << std::endl;
		// // auto daily_std = Daily_Modulation_Standard(points, DM_standard, SHM, graphene, 1e5);
		// // libphysica::Export_Table("Daily_Modulation_Standard.txt", daily_std, {1.0, dim});
		// std::cout << "Tabulate daily modulation for nreft" << std::endl;
		// auto daily_nreft = Daily_Modulation_NREFT(points, DM_nreft, SHM, graphene, 2e6);
		// libphysica::Export_Table("Daily_Modulation_NREFT.txt", daily_nreft, {1.0, dim});
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}