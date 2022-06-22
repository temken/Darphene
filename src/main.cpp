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

#include "DM_Particle_NREFT.hpp"
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

	// Initialize DM halo model (values from 2105.00599)
	double rho_DM  = 0.3 * GeV / cm / cm / cm;
	double v0	   = 238.0 * km / sec;
	double v_earth = 250.552 * km / sec;
	double v_esc   = 544.0 * km / sec;
	obscura::Standard_Halo_Model SHM(rho_DM, v0, v_earth, v_esc);
	double t = 0.0 * hr;
	SHM.Set_Observer_Velocity(Earth_Velocity(t, v_earth));
	SHM.Print_Summary();

	// Initialize graphene
	double work_function = 5 * eV;
	Graphene graphene(work_function);

	obscura::DM_Particle_SI DM_standard(100 * MeV);
	DM_standard.Set_Sigma_Electron(1.0e-36 * cm * cm);
	DM_standard.Print_Summary();

	DM_Particle_NREFT DM_nreft(100.0 * MeV, 0.5);
	// DM.Set_Coupling(3, 1.0, "Contact");
	// DM.Set_Coupling(5, 0.1, "Power", 3);
	// DM.Set_Coupling(7, 0.1, "General", 100 * MeV);
	DM_nreft.Set_Cross_Section(1, pb, "Contact");
	DM_nreft.Print_Summary();

	int band				 = 0;
	double R_standard		 = R_Total_Standard(DM_standard, SHM, graphene, band, 1e5);			 //, band);
	double R_standard_simple = R_Total_Standard_Simple(DM_standard, SHM, graphene, band, 1e5);	 //, band);
	double R_nreft			 = R_Total_NREFT(DM_nreft, SHM, graphene, 0, 1e6);					 //, band);

	std::cout << In_Units(R_standard, 1.0 / kg / year) << " " << In_Units(R_standard_simple, 1.0 / kg / year) << " " << In_Units(R_nreft, 1.0 / kg / year) << std::endl;

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}