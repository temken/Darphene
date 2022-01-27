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

	std::vector<std::vector<std::vector<double>>> energy_bands;
	std::vector<std::vector<std::vector<double>>> spectra;
	for(int i = 0; i < 50; i++)
	{
		std::cout << i << std::endl;
		Graphene graphene;
		auto energy_band = graphene.Energy_Bands(100);
		energy_bands.push_back(energy_band);

		std::vector<double> energies = libphysica::Log_Space(2.0e-1 * eV, 250 * eV, 50);
		std::vector<std::vector<double>> spectrum;
		for(auto& E_e : energies)
		{
			double value = 0.0;
			for(int band = 0; band < 4; band++)
				value += dR_dlnE_corrected(E_e, DM, SHM, graphene, band);
			spectrum.push_back({E_e, value});
		}
		spectra.push_back(spectrum);
	}
	std::cout << energy_bands.size() << std::endl;
	// std::cout << spectra.size() << std::endl;
	std::cout << energy_bands[0].size() << std::endl;
	// std::cout << spectra[0].size() << std::endl;

	// Compute standard deviations of energy bands and spectra
	// 1. Energy bands
	std::ofstream f;
	f.open("Variation_Energy_Bands_0_1.txt");
	for(int i = 0; i < energy_bands[0].size(); i++)
	{
		std::vector<double> values_pi;
		std::vector<double> values_s1;
		std::vector<double> values_s2;
		std::vector<double> values_s3;
		for(auto& band : energy_bands)
		{
			values_pi.push_back(band[i][1]);
			values_s1.push_back(band[i][3]);
			values_s2.push_back(band[i][4]);
			values_s3.push_back(band[i][5]);
		}
		double mu_pi	= libphysica::Arithmetic_Mean(values_pi);
		double sigma_pi = libphysica::Standard_Deviation(values_pi);
		double mu_s1	= libphysica::Arithmetic_Mean(values_s1);
		double sigma_s1 = libphysica::Standard_Deviation(values_s1);
		double mu_s2	= libphysica::Arithmetic_Mean(values_s2);
		double sigma_s2 = libphysica::Standard_Deviation(values_s2);
		double mu_s3	= libphysica::Arithmetic_Mean(values_s3);
		double sigma_s3 = libphysica::Standard_Deviation(values_s3);
		f << energy_bands[0][i][0] / keV << "\t" << mu_pi / eV << "\t" << sigma_pi / eV
		  << "\t" << mu_s1 / eV << "\t" << sigma_s1 / eV
		  << "\t" << mu_s2 / eV << "\t" << sigma_s2 / eV
		  << "\t" << mu_s3 / eV << "\t" << sigma_s3 / eV
		  << std::endl;
	}
	f.close();

	// Spectra
	f.open("Variation_Spectrum_0_1.txt");
	for(int i = 0; i < spectra[0].size(); i++)
	{
		std::vector<double> values;
		for(auto& spectrum : spectra)
			values.push_back(spectrum[i][1]);
		double mu	 = libphysica::Arithmetic_Mean(values);
		double sigma = libphysica::Standard_Deviation(values);
		f << In_Units(spectra[0][i][0], eV) << "\t" << In_Units(mu, 1.0 / kg / year) << "\t" << In_Units(sigma, 1.0 / kg / year) << std::endl;
	}
	f.close();
	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	return 0;
}