#ifndef __DM_Particle_NREFT_hpp_
#define __DM_Particle_NREFT_hpp_

#include <Eigen/Geometry>

#include "obscura/DM_Particle.hpp"

namespace graphene
{

// 1. Small class for different kind of DM form factors that encapsulate additional q-dependences in the amplitude
class DM_Form_Factor
{
  private:
	std::string form_factor_type;
	double m_Mediator;
	int q_power;
	double qRef;

  public:
	DM_Form_Factor();
	DM_Form_Factor(std::string& type, double parameter = 0.0);

	void Print_Summary(int rank = 0) const;

	double operator()(double q);
};

// 2. Class for the DM particle interacting with electrons via effective operators O_i
class DM_Particle_NREFT : public obscura::DM_Particle
{
  protected:
	std::vector<double> couplings;
	std::vector<DM_Form_Factor> DM_form_factors;

  public:
	DM_Particle_NREFT();
	DM_Particle_NREFT(double m, double s = 0.5);

	virtual double Get_Interaction_Parameter(std::string target) const override;
	virtual void Set_Interaction_Parameter(double par, std::string target) override;

	void Set_Coupling(int index, double value, std::string form_factor = "Contact", double param = 0.0);
	void Reset_All_Couplings();

	double Squared_Amplitude_Electron(const Eigen::Vector3d& qVec, const Eigen::Vector3d& velDM, const Eigen::Vector3d& kPrime);

	virtual double Sigma_Total_Electron(double vDM, double param = -1.0) override;

	virtual void Print_Summary(int rank = 0) const override;
};

}	// namespace graphene
#endif