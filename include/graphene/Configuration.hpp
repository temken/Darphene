#ifndef __Configuration_hpp_
#define __Configuration_hpp_

#include "obscura/Configuration.hpp"

#include "graphene/Graphene.hpp"

namespace graphene
{
// 1. Configuration class for input file, which extends the obscura::Configuration class.

class Configuration : public obscura::Configuration
{
  protected:
	virtual void Construct_DM_Particle() override;
	void Construct_DM_Particle_Dark_Photon(double mDM);
	void Construct_DM_Particle_NREFT(double mDM);
	void Construct_DM_Particle_Electric_Dipole(double mDM);
	void Construct_DM_Particle_Magnetic_Dipole(double mDM);
	void Construct_DM_Particle_Anapole(double mDM);

	void Import_Graphene_Parameters();

	std::vector<double> NREFT_couplings;
	std::vector<std::string> NREFT_form_factors;
	std::vector<double> NREFT_parameters;

  public:
	Graphene graphene;
	std::string run_modus;
	unsigned int MC_points, grid_points;

	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int mpi_rank = 0) override;
};

}	// namespace graphene

#endif