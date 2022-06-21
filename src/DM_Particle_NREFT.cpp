#include "DM_Particle_NREFT.hpp"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

namespace graphene
{

using namespace libphysica::natural_units;

// 1. Small class for different kind of DM form factors that encapsulate additional q-dependences in the amplitude

DM_Form_Factor::DM_Form_Factor()
: form_factor_type("Contact"), m_Mediator(0.0), q_power(0), qRef(aEM * mElectron)
{
}

DM_Form_Factor::DM_Form_Factor(std::string& type, double parameter)
: form_factor_type(type), qRef(aEM * mElectron)
{
	if(form_factor_type == "General")
		m_Mediator = parameter;
	else if(form_factor_type == "Power")
		q_power = parameter;
}

void DM_Form_Factor::Print_Summary(int rank) const
{
	if(rank == 0)
	{
		std::cout << form_factor_type;
		if(form_factor_type == "General")
			std::cout << "(" << libphysica::Round(m_Mediator / MeV) << " MeV)";
		else if(form_factor_type == "Power")
			std::cout << "(" << q_power << ")";
	}
}

double DM_Form_Factor::operator()(double q)
{
	if(form_factor_type == "Contact")
		return 1.0;
	else if(form_factor_type == "Long-Range")
		return qRef * qRef / q / q;
	else if(form_factor_type == "General")
		return (qRef * qRef + m_Mediator * m_Mediator) / (q * q + m_Mediator * m_Mediator);
	else if(form_factor_type == "Power")
		return std::pow(q / qRef, q_power);
	else
	{
		std::cerr << "Error in DM_Form_Factor::operator(): Form factor " << form_factor_type << " not implemented." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// 2. Class for the DM particle interacting with electrons via effective operators O_i

DM_Particle_NREFT::DM_Particle_NREFT()
: obscura::DM_Particle(100.0 * MeV, 0.5), couplings(std::vector<double>(15, 0.0)), DM_form_factors(std::vector<DM_Form_Factor>(15))
{
}

DM_Particle_NREFT::DM_Particle_NREFT(double m, double s)
: obscura::DM_Particle(m, s), couplings(std::vector<double>(15, 0.0)), DM_form_factors(std::vector<DM_Form_Factor>(15))
{
}

double DM_Particle_NREFT::Get_Interaction_Parameter(std::string target) const
{
	return 0.0;
}
void DM_Particle_NREFT::Set_Interaction_Parameter(double par, std::string target)
{
}

void DM_Particle_NREFT::Set_Coupling(int index, double value, std::string form_factor, double param)
{
	if(index < 1 || index > 15 || index == 2)
	{
		std::cerr << "Error in DM_Particle_NREFT::Set_Coupling: Operator " << index << " is undefined." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	couplings[index - 1]	   = value;
	DM_form_factors[index - 1] = DM_Form_Factor(form_factor, param);
}

void DM_Particle_NREFT::Reset_All_Couplings()
{
	couplings		= std::vector<double>(15, 0.0);
	DM_form_factors = std::vector<DM_Form_Factor>(15);
}

double DM_Particle_NREFT::Squared_Amplitude_Electron(const Eigen::Vector3d& qVec, const Eigen::Vector3d& velDM, const Eigen::Vector3d& kPrime)
{
	// Kinematic quantities

	// Effective couplings including form factors

	// Squared amplitude
	double M2 = 0.0;
	return M2;
}

double DM_Particle_NREFT::Sigma_Total_Electron(double vDM, double param)
{
	return 0.0;
}

void DM_Particle_NREFT::Print_Summary(int rank) const
{
	if(rank == 0)
	{
		Print_Summary_Base(rank);
		std::cout << "\nOperator\tValue\tForm Factor" << std::endl
				  << "----------------------------------------" << std::endl;
		for(int op = 1; op < 16; op++)
		{
			if(op == 2)
				continue;
			else if(couplings[op - 1] != 0.0)
			{
				std::cout << op << "\t\t" << libphysica::Round(couplings[op - 1]) << "\t";
				DM_form_factors[op - 1].Print_Summary(rank);
				std::cout << std::endl;
			}
		}
		std::cout << "----------------------------------------" << std::endl;
	}
}

}	// namespace graphene