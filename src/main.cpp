#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Numerics.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"
#include "DM_Particle_Standard.hpp"

#include "Direct_Detection_Graphene.hpp"
#include "Graphene.hpp"
#include "Hydrogenic_Wavefunctions.hpp"
#include "Multidimensional_Integration.hpp"
#include "Vegas.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace graphene;
// using namespace boost::math::quadrature;

int main(int argc, char* argv[])
{
	//Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << std::endl;
	////////////////////////////////////////////////////////////////////////

	obscura::Standard_Halo_Model SHM(0.4 * GeV / cm / cm / cm, 220.0 * km / sec, 232.0 * km / sec, 550.0 * km / sec);
	obscura::DM_Particle_SI DM(100 * MeV);
	Graphene graphene;

	// std::vector<double> energies = libphysica::Log_Space(0.2 * eV, 250 * eV, 50);
	// for(int i = 0; i < 1; i++)
	// {
	// 	std::ofstream f;
	// 	f.open("Graphene_Spectrum_" + std::to_string(i) + ".txt");
	// 	for(auto& E : energies)
	// 		f << E / eV << "\t" << kg * year * 2.0 * perform_integral_vegas(E, DM, SHM, graphene, i) << std::endl;
	// 	f.close();
	// }

	double E_e = 100 * eV;
	// for(int i = 0; i < 4; i++)
	// 	std::cout << kg * year * perform_integral(E_e, DM, SHM, graphene, i) << std::endl;
	for(int i = 0; i < 1; i++)
	{
		std::cout << perform_integral(E_e, DM, SHM, graphene, i) << std::endl;
		std::cout << 2.0 * perform_integral_vegas(E_e, DM, SHM, graphene, i) << std::endl
				  << std::endl;
	}

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}