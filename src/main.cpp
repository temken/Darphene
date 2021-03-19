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

#include "Graphene.hpp"
#include "Hydrogenic_Wavefunctions.hpp"
#include "Multidimensional_Integration.hpp"
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

	double aCC = 1.42 * Angstrom;
	double a   = aCC * sqrt(3.0);
	Eigen::Vector3d high_symmetry_point_G({0.0, 0.0, 0.0});
	Eigen::Vector3d high_symmetry_point_M({2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0});
	Eigen::Vector3d high_symmetry_point_K({2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0});

	Graphene graphene;

	Eigen::Vector3d rVec({Bohr_Radius, 2 * Bohr_Radius, 3 * Bohr_Radius});
	Eigen::Vector3d lVec({1.0 / Bohr_Radius, 2 / Bohr_Radius, 3.0 / Bohr_Radius});
	// for(auto entry : graphene.Energy_Dispersion_Sigma(high_symmetry_point_G))
	// 	std::cout << entry / eV << std::endl;

	std::complex<double> psi		  = graphene.Wavefunction_Pi(rVec, lVec);
	std::complex<double> psi_analytic = graphene.Wavefunction_Pi_Analytic(rVec, lVec);
	std::cout << psi << "\t" << psi_analytic << std::endl;
	std::cout << std::abs(psi) << "\t" << std::abs(psi_analytic) << std::endl;

	double energy		   = graphene.Energy_Dispersion_Pi(high_symmetry_point_M)[1];
	double energy_analytic = graphene.Energy_Dispersion_Pi_Analytic(high_symmetry_point_M)[1];
	std::cout << energy / eV << "\t" << energy_analytic / eV << std::endl;

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}