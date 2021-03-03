#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include <Eigen/Eigenvalues>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

#include "version.hpp"

using namespace libphysica::natural_units;

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

	Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> ges;
	Eigen::MatrixXf A = Eigen::MatrixXf::Random(4, 4);
	Eigen::MatrixXf B = Eigen::MatrixXf::Random(4, 4);
	ges.compute(A, B);
	std::cout << "The (complex) numerators of the generalzied eigenvalues are: " << ges.alphas().transpose() << std::endl;
	std::cout << "The (real) denominatore of the generalzied eigenvalues are: " << ges.betas().transpose() << std::endl;
	std::cout << "The (complex) generalzied eigenvalues are (alphas./beta): " << ges.eigenvalues().transpose() << std::endl;

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}