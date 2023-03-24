#include "graphene/DM_Particle_NREFT.hpp"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

namespace Darphene
{

using namespace libphysica::natural_units;

// 1. Small class for different kind of DM form factors that encapsulate additional q-dependences in the amplitude

DM_Form_Factor::DM_Form_Factor()
: form_factor_type("Contact"), m_Mediator(0.0), q_power(0), qRef(aEM * mElectron)
{
}

DM_Form_Factor::DM_Form_Factor(const std::string& type, double parameter)
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
	if(form_factor_type == "Contact" || form_factor_type == "C")
		return 1.0;
	else if(form_factor_type == "Long-Range" || form_factor_type == "L")
		return qRef * qRef / q / q;
	else if(form_factor_type == "General")
		return (qRef * qRef + m_Mediator * m_Mediator) / (q * q + m_Mediator * m_Mediator);
	else if(form_factor_type == "Power")
		return std::pow(q / qRef, q_power);
	else
	{
		std::cerr << "\033[1;31mError\033[0m in DM_Form_Factor::operator(): Form factor " << form_factor_type << " not implemented." << std::endl;
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

void DM_Particle_NREFT::Set_Coupling(int index, double value, const std::string& form_factor, double param)
{
	if(index < 1 || index > 15 || index == 2)
	{
		std::cerr << "\033[1;31mError\033[0m in DM_Particle_NREFT::Set_Coupling: Operator " << index << " is undefined." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	couplings[index - 1]	   = value;
	DM_form_factors[index - 1] = DM_Form_Factor(form_factor, param);
}

void DM_Particle_NREFT::Set_Cross_Section(int index, double sigma, const std::string& form_factor, double param)
{
	double c = 4.0 * std::sqrt(M_PI * sigma) * mass * mElectron / libphysica::Reduced_Mass(mass, mElectron);
	Set_Coupling(index, c, form_factor, param);
}

void DM_Particle_NREFT::Reset_All_Couplings()
{
	couplings		= std::vector<double>(15, 0.0);
	DM_form_factors = std::vector<DM_Form_Factor>(15);
}

double DM_Particle_NREFT::Response_Function(const Eigen::Vector3d& qVec, const Eigen::Vector3d& velDM, const Eigen::Vector3d& kPrime)
{
	// Kinematic quantities
	double q				  = qVec.norm();
	Eigen::Vector3d q_minus_k = (qVec - kPrime) / mElectron;
	Eigen::Vector3d vel_perp  = velDM - qVec / 2.0 / libphysica::Reduced_Mass(mass, mElectron);
	double v_perp_sq		  = vel_perp.squaredNorm();
	double v_perp_dot_q		  = vel_perp.dot(qVec) / mElectron;
	double q_me_sq			  = q * q / mElectron / mElectron;
	double v_perp_cross_q_sq  = q_me_sq * v_perp_sq - v_perp_dot_q * v_perp_dot_q;
	double q_minus_k_sq		  = q_minus_k.squaredNorm();
	double q_qk_prime		  = q_minus_k.dot(qVec) / mElectron;

	// Effective couplings including form factors
	std::vector<double> C(15, 0.0);
	for(int i = 0; i < 15; i++)
		C[i] = couplings[i] * DM_form_factors[i](q);

	// Squared amplitude
	double M2_1 = C[0] * C[0] + C[2] * C[2] / 4. * v_perp_cross_q_sq + C[6] * C[6] / 4. * v_perp_sq + C[9] * C[9] / 4. * q_me_sq;
	double M2_2 = 3. * C[3] * C[3] + (4. * C[4] * C[4] - 2. * C[11] * C[14]) * v_perp_cross_q_sq + C[5] * C[5] * q_me_sq * q_me_sq + (4. * C[7] * C[7] + 2. * C[11] * C[11]) * v_perp_sq + (2. * C[8] * C[8] + 4. * C[10] * C[10] + 2. * C[3] * C[5]) * q_me_sq + (C[12] * C[12] + C[13] * C[13]) * q_me_sq * v_perp_sq + C[14] * C[14] * q_me_sq * v_perp_cross_q_sq + 2. * C[12] * C[13] * v_perp_dot_q * v_perp_dot_q;
	double M2	= M2_1 + spin * (spin + 1.0) / 12.0 * M2_2;

	double first_bracket  = (C[2] * C[2] / 2. * (v_perp_dot_q * qVec / mElectron - q_me_sq * vel_perp) - C[6] * C[6] / 2. * vel_perp).dot(q_minus_k);
	double second_bracket = ((4. * C[4] * C[4] + C[14] * C[14] * q_me_sq) * (v_perp_dot_q * qVec / mElectron - q_me_sq * vel_perp) - (4. * C[7] * C[7] + 2. * C[11] * C[11] + (C[12] * C[12] + C[13] * C[13]) * q_me_sq) * vel_perp - 2. * C[11] * C[14] * (v_perp_dot_q * qVec / mElectron - q_me_sq * vel_perp) - 2. * C[12] * C[13] * v_perp_dot_q * qVec / mElectron).dot(q_minus_k);
	double me_real_part	  = first_bracket + spin * (spin + 1.0) / 6.0 * second_bracket;

	double first_term	   = (C[2] * C[2] / 4. * q_me_sq + C[6] * C[6] / 4.) * q_minus_k_sq - C[2] * C[2] / 4. * q_qk_prime * q_qk_prime;
	double bracket		   = q_minus_k_sq * ((4. * C[4] * C[4] + C[12] * C[12] + C[13] * C[13] - 2. * C[11] * C[14]) * q_me_sq + 4. * C[7] * C[7] + 2. * C[11] * C[11] + C[14] * C[14] * q_me_sq * q_me_sq) + q_qk_prime * q_qk_prime * (-4. * C[4] * C[4] - C[14] * C[14] * q_me_sq + 2. * C[11] * C[14] + 2. * C[12] * C[13]);
	double me_squared_part = first_term + spin * (spin + 1.0) / 12.0 * bracket;

	return M2 + me_real_part + me_squared_part;
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

DM_Particle_NREFT DM_Dark_Photon(double mDM, double sigma_e, const std::string& form_factor, double mMediator)
{
	DM_Particle_NREFT DM(mDM);
	DM.Set_Cross_Section(1, sigma_e, form_factor, mMediator);
	return DM;
}

DM_Particle_NREFT DM_Electric_Dipole(double mDM, double g_over_lambda)
{
	DM_Particle_NREFT DM(mDM);
	double qRef = aEM * mElectron;
	double c11	= 16.0 * Elementary_Charge * mDM * mElectron * mElectron * g_over_lambda / qRef / qRef;
	DM.Set_Coupling(11, c11, "Long-Range");
	return DM;
}

DM_Particle_NREFT DM_Magnetic_Dipole(double mDM, double g_over_lambda)
{
	DM_Particle_NREFT DM(mDM);
	double qRef = aEM * mElectron;
	double c1	= 4.0 * Elementary_Charge * mElectron * g_over_lambda;
	double c4	= 4.0 * Elementary_Charge * mDM * g_over_lambda;
	double c5	= 16.0 * Elementary_Charge * mDM * mElectron * mElectron * g_over_lambda / qRef / qRef;
	double c6	= -c5;
	DM.Set_Coupling(1, c1, "Contact");
	DM.Set_Coupling(4, c4, "Contact");
	DM.Set_Coupling(5, c5, "Long-Range");
	DM.Set_Coupling(6, c6, "Long-Range");
	return DM;
}

DM_Particle_NREFT DM_Anapole(double mDM, double g_over_lambda_2)
{
	DM_Particle_NREFT DM(mDM);

	double c8 = 8.0 * Elementary_Charge * mElectron * mDM * g_over_lambda_2;
	double c9 = -c8;
	DM.Set_Coupling(8, c8, "Contact");
	DM.Set_Coupling(9, c9, "Contact");
	return DM;
}

}	// namespace Darphene