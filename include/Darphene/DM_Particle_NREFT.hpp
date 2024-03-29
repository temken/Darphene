#ifndef __DM_Particle_NREFT_hpp_
#define __DM_Particle_NREFT_hpp_

#include <Eigen/Geometry>

#include "obscura/DM_Particle.hpp"

namespace Darphene
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
	DM_Form_Factor(const std::string& type, double parameter = 0.0);

	void Print_Summary(int rank = 0) const;

	double operator()(double q) const;
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

	void Set_Coupling(int index, double value, const std::string& form_factor = "Contact", double param = 0.0);
	void Set_Cross_Section(int index, double sigma, const std::string& form_factor = "Contact", double param = 0.0);
	void Reset_All_Couplings();

	virtual double Sigma_Electron() const override;

	double Response_Function(const Eigen::Vector3d& qVec, const Eigen::Vector3d& velDM, const Eigen::Vector3d& kPrime) const;

	virtual void Print_Summary(int rank = 0) const override;
};

// 3. Class for benchmark models
extern DM_Particle_NREFT DM_Dark_Photon(double mDM, double sigma_e, const std::string& form_factor = "Contact", double mMediator = 0.0);
extern DM_Particle_NREFT DM_Electric_Dipole(double mDM, double g_over_lambda);
extern DM_Particle_NREFT DM_Magnetic_Dipole(double mDM, double g_over_lambda);
extern DM_Particle_NREFT DM_Anapole(double mDM, double g_over_lambda_2);

}	// namespace Darphene

#endif