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

	obscura::Standard_Halo_Model SHM(0.4 * GeV / cm / cm / cm, 220.0 * km / sec, 232.0 * km / sec, 550.0 * km / sec);
	obscura::DM_Particle_SI DM(100 * MeV);
	DM.Set_Sigma_Electron(1.0e-37 * cm * cm);
	Graphene graphene;

	double E_e = 10 * eV;
	int band   = 0;
	std::cout << E_e / eV << "\t" << In_Units(dR_dlnE(E_e, DM, SHM, graphene, band), 1.0 / kg / year) << std::endl;
	std::cout << E_e / eV << "\t" << In_Units(dR_dlnE_corrected(E_e, DM, SHM, graphene, band), 1.0 / kg / year) << std::endl;
	// std::cout << E_e / eV << "\t" << In_Units(dR_dlnE_corrected(E_e, DM, SHM, graphene, band), 1.0 / kg / year) << std::endl;

	// Eigen::Vector3d rVec(Bohr_Radius, 2.0 * Bohr_Radius, 3.0 * Bohr_Radius);
	// Eigen::Vector3d kVec(1.0 / Bohr_Radius, 2.0 / Bohr_Radius, 3.0 / Bohr_Radius);
	// Eigen::Vector3d lVec(1.0 / Bohr_Radius, 2.0 / Bohr_Radius, 0.0);

	// std::cout << graphene.Wavefunction_Pi(rVec, lVec) << std::endl;
	// std::cout << graphene.Wavefunction_Pi_Analytic(rVec, lVec) << std::endl
	// 		  << std::endl;

	// std::cout << graphene.Wavefunction_Momentum_Pi(kVec, lVec) << std::endl;
	// std::cout << graphene.Wavefunction_Momentum_Pi_Analytic(kVec, lVec) << std::endl;

	// std::cout << DM.mass / MeV << std::endl;

	std::vector<double> energies = libphysica::Log_Space(2.0e-1 * eV, 200 * eV, 100);
	std::ofstream f;
	f.open("Corrected_Spectrum_Pi_100_MeV.txt");
	for(auto& E_e : energies)
	{
		std::cout << E_e / eV << std::endl;
		f << In_Units(E_e, eV) << "\t" << In_Units(dR_dlnE_corrected(E_e, DM, SHM, graphene, 0), 1.0 / kg / year) << std::endl;
	}
	f.close();

	// std::vector<double> energies = libphysica::Log_Space(2.0e-1 * eV, 30 * eV, 50);
	// std::ofstream f;
	// f.open("Spectrum_Emken_10_MeV.txt");
	// for(auto& E_e : energies)
	// {
	// 	std::cout << E_e / eV << std::endl;
	// 	f << In_Units(E_e, eV);
	// 	for(int band = 0; band < 4; band++)
	// 		f << "\t" << In_Units(dR_dlnE(E_e, DM, SHM, graphene, band), 1.0 / kg / year);
	// 	f << std::endl;
	// }
	// f.close();

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}