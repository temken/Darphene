#include "graphene/Configuration.hpp"

#include <libconfig.h++>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"

#include "graphene/DM_Particle_NREFT.hpp"
#include "graphene/Direct_Detection_Standard.hpp"
#include "version.hpp"

namespace DarPhene
{
using namespace libconfig;
using namespace libphysica::natural_units;

void Configuration::Construct_DM_Particle()
{
	double DM_mass, DM_spin, DM_fraction;
	// 3.1 General properties
	try
	{
		DM_mass = config.lookup("DM_mass");
		DM_mass *= MeV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_mass' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_spin = config.lookup("DM_spin");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_spin' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	try
	{
		DM_fraction = config.lookup("DM_fraction");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_fraction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 3.2 DM interactions
	std::string DM_interaction;
	try
	{
		DM_interaction = config.lookup("DM_interaction").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'DM_interaction' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// DM Models
	if(DM_interaction == "Dark-Photon")
		Configuration::Construct_DM_Particle_Dark_Photon(DM_mass);
	else if(DM_interaction == "NREFT")
		Configuration::Construct_DM_Particle_NREFT(DM_mass);
	else if(DM_interaction == "Electric-Dipole")
		Configuration::Construct_DM_Particle_Electric_Dipole(DM_mass);
	else if(DM_interaction == "Magnetic-Dipole")
		Configuration::Construct_DM_Particle_Magnetic_Dipole(DM_mass);
	else if(DM_interaction == "Anapole")
		Configuration::Construct_DM_Particle_Anapole(DM_mass);
	else
	{
		std::cerr << "\033[1;31mError\033[0m in graphene::Configuration::Construct_DM_Particle(): 'DM_interaction' setting " << DM_interaction << " in configuration file not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM_NREFT->Set_Spin(DM_spin);
	DM_NREFT->Set_Fractional_Density(DM_fraction);
	DM_NREFT->Set_Low_Mass_Mode(true);
}

void Configuration::Construct_DM_Particle_Dark_Photon(double mDM)
{
	// DM form factor
	std::string DM_form_factor;
	double DM_mediator_mass = -1.0;
	try
	{
		DM_form_factor = config.lookup("DM_form_factor").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "\033[1;31mError\033[0m in DaMaSCUS_SUN::Configuration::Construct_DM_Particle_Dark_Photon(): No 'DM_form_factor' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(DM_form_factor == "General")
	{
		try
		{
			DM_mediator_mass = config.lookup("DM_mediator_mass");
			DM_mediator_mass *= MeV;
		}
		catch(const SettingNotFoundException& nfex)
		{
			std::cerr << "\033[1;31mError\033[0m in Configuration::Construct_DM_Particle_Dark_Photon(): No 'DM_mediator_mass' setting in configuration file." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	double DM_cross_section_electron;
	try
	{
		DM_cross_section_electron = config.lookup("DM_cross_section_electron");
		DM_cross_section_electron *= cm * cm;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "\033[1;31mError\033[0m in Configuration::Construct_DM_Particle_Dark_Photon(): No 'DM_cross_section_electron' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM_NREFT = new DM_Particle_NREFT(DM_Dark_Photon(mDM, DM_cross_section_electron, DM_form_factor, DM_mediator_mass));
}

void Configuration::Construct_DM_Particle_NREFT(double mDM)
{
	DM_NREFT	= new DM_Particle_NREFT(mDM);
	int entries = config.lookup("NREFT_couplings").getLength();
	for(int i = 0; i < entries; i++)
	{
		int op					= config.lookup("NREFT_couplings")[i][0];
		double coupling			= config.lookup("NREFT_couplings")[i][1];
		std::string form_factor = config.lookup("NREFT_couplings")[i][2].c_str();
		double param			= config.lookup("NREFT_couplings")[i][3];
		if(form_factor == "General")
			param *= MeV;
		dynamic_cast<DM_Particle_NREFT*>(DM_NREFT)->Set_Coupling(op, coupling, form_factor, param);
	}
}

void Configuration::Construct_DM_Particle_Electric_Dipole(double mDM)
{
	double DM_coupling;
	try
	{
		DM_coupling = config.lookup("DM_coupling");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "\033[1;31mError\033[0m in Configuration::Construct_DM_Particle_Electric_Dipole(): No 'DM_coupling' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM_NREFT = new DM_Particle_NREFT(DM_Electric_Dipole(mDM, DM_coupling));
}

void Configuration::Construct_DM_Particle_Magnetic_Dipole(double mDM)
{
	double DM_coupling;
	try
	{
		DM_coupling = config.lookup("DM_coupling");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "\033[1;31mError\033[0m in Configuration::Construct_DM_Particle_Magnetic_Dipole(): No 'DM_coupling' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM_NREFT = new DM_Particle_NREFT(DM_Magnetic_Dipole(mDM, DM_coupling));
}

void Configuration::Construct_DM_Particle_Anapole(double mDM)
{
	double DM_coupling;
	try
	{
		DM_coupling = config.lookup("DM_coupling");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "\033[1;31mError\033[0m in Configuration::Construct_DM_Particle_Anapole(): No 'DM_coupling' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DM_NREFT = new DM_Particle_NREFT(DM_Anapole(mDM, DM_coupling));
}

void Configuration::Import_Graphene_Parameters()
{
	try
	{
		run_modus = config.lookup("run_modus").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'run_modus' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		MC_points = config.lookup("MC_points");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'MC_points' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		grid_points = config.lookup("grid_points");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'grid_points' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Graphene - Wavefunctions
	try
	{
		carbon_wavefunctions = config.lookup("carbon_wave_functions").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'carbon_wave_functions' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Graphene - Work function
	try
	{
		graphene_work_function = config.lookup("work_function");
		graphene_work_function *= eV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'work_function' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Time of the day (t=0 is DM wind from top, t=12 DM wind parallel to graphene)
	try
	{
		time = config.lookup("time");
		time *= hr;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'time' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
{
	cfg_file	 = cfg_filename;
	results_path = "./";

	// 1. Read the cfg file.
	Read_Config_File();

	// 2. Find the run ID, create a folder and copy the cfg file.
	Initialize_Result_Folder(MPI_rank);

	// 3. DM particle
	Construct_DM_Particle();

	// 4. DM Distribution
	Construct_DM_Distribution();

	// 5. Computation of exclusion limits
	Initialize_Parameters();

	// 6. Graphene specific parameters
	Import_Graphene_Parameters();

	// Rotate the vel_Earth vector.
	double vEarth = dynamic_cast<obscura::Standard_Halo_Model*>(DM_distr)->Get_Observer_Velocity().Norm();
	dynamic_cast<obscura::Standard_Halo_Model*>(DM_distr)->Set_Observer_Velocity(Earth_Velocity(time, vEarth));
}

void Configuration::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << "Summary of graphene configuration" << std::endl
				  << std::endl
				  << "Config file:\t" << cfg_file << std::endl
				  << "ID:\t\t" << ID << std::endl;
		DM_NREFT->Print_Summary(mpi_rank);
		DM_distr->Print_Summary(mpi_rank);
		std::cout << "Graphene options" << std::endl
				  << "\tRun modus:\t\t\t" << run_modus << std::endl
				  << "\tMC points:\t\t\t" << MC_points << std::endl
				  << "\tGrid points:\t\t\t" << grid_points << std::endl
				  << "\tWave functions:\t\t\t" << carbon_wavefunctions << std::endl
				  << "\tGraphene work function [eV]:\t" << graphene_work_function / eV << std::endl
				  << "\tTime of day [hr]:\t\t" << time / hr << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DarPhene