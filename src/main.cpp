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

#include "Direct_Detection_Graphene.hpp"
#include "Graphene.hpp"
#include "Hydrogenic_Wavefunctions.hpp"
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

	// Initialize DM halo model
	double rho_DM  = 0.4 * GeV / cm / cm / cm;
	double v0	   = 220.0 * km / sec;
	double v_earth = 232.0 * km / sec;
	double v_esc   = 544.0 * km / sec;
	obscura::Standard_Halo_Model SHM(rho_DM, v0, v_earth, v_esc);
	double t = 0.2 * day;
	SHM.Set_Observer_Velocity(Earth_Velocity(t, v_earth));
	SHM.Print_Summary();

	// Initialize graphene
	double work_function = 5 * eV;
	Graphene graphene(work_function);

	// Initialize DM particle
	obscura::DM_Particle_SI DM(400 * MeV);
	DM.Set_Sigma_Electron(1.0e-37 * cm * cm);

	// Compute recoil spectrum
	int points					  = 40;
	std::string velocity_integral = "Full";
	auto spectrum				  = Tabulate_dR_dlnE(points, DM, SHM, graphene, velocity_integral);
	libphysica::Export_Table("Spectrum_mDM_400MeV_Simplified.txt", spectrum, {eV, 1.0 / kg / year, 1.0 / kg / year, 1.0 / kg / year, 1.0 / kg / year, 1.0 / kg / year});

	// Total rate over the course of a day
	double R_simple = R_Total_simplified(DM, SHM, graphene, 0);
	std::cout << "\tR_simple = " << In_Units(R_simple, 1.0 / kg / year) << std::endl;
	auto t_list = libphysica::Linear_Space(0.0, 1.0, 20);
	for(auto& t : t_list)
	{
		SHM.Set_Observer_Velocity(Earth_Velocity(t * day, v_earth));
		double R = R_Total(DM, SHM, graphene, 0);
		std::cout << "t = " << t << "\tR = " << In_Units(R, 1.0 / kg / year) << std::endl;
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}