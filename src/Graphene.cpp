#include "Graphene.hpp"

#include <algorithm>

#include "libphysica/Utilities.hpp"

#include "Hydrogenic_Wavefunctions.hpp"

namespace graphene
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

std::complex<double> Graphene::f_aux(const Eigen::Vector3d& kVec) const
{
	return exp(1i * kVec[0] * a / sqrt(3.0)) + 2.0 * exp(-1i * kVec[0] * a / 2.0 / sqrt(3.0)) * cos(a * kVec[1] / 2.0);
}

Eigen::Vector3d Graphene::Path_1BZ(double k) const
{
	if(k >= 0.0 && k < 4.0 * M_PI / 3.0 / a)
		return {2.0 * M_PI / sqrt(3.0) / a - sqrt(3.0) / 2.0 * k, -k / 2.0 + 2.0 * M_PI / 3.0 / a, 0.0};
	else if(k >= 4.0 * M_PI / 3.0 / a && k < 2.0 * M_PI / 3.0 / a * (2.0 + sqrt(3.0)))
		return {k - 4.0 / 3.0 * M_PI / a, 0.0, 0.0};
	else if(k >= 2.0 * M_PI / 3.0 / a * (2.0 + sqrt(3.0)) && k <= 2.0 * M_PI / 3.0 / a * (3.0 + sqrt(3.0)))
		return {2.0 * M_PI / sqrt(3.0) / a, k - 2.0 / 3.0 * M_PI * (2.0 + sqrt(3.0)) / a, 0.0};
	else
	{
		std::cerr << "Error in Graphene::Path_1BZ(): Momentum k = " << k << " out of bound." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

Eigen::MatrixXcd Graphene::S_Matrix_Pi(const Eigen::Vector3d& kVec) const
{
	Eigen::MatrixXcd S(2, 2);
	S << 1.0, s * f_aux(kVec), s * std::conj(f_aux(kVec)), 1.0;
	return S;
}

Eigen::MatrixXcd Graphene::H_Matrix_Pi(const Eigen::Vector3d& kVec) const
{
	Eigen::MatrixXcd H(2, 2);
	H << epsilon_2p, t * f_aux(kVec), t * std::conj(f_aux(kVec)), epsilon_2p;
	return H;
}

Eigen::MatrixXcd Graphene::S_Matrix_Sigma(const Eigen::Vector3d& kVec) const
{
	Eigen::MatrixXcd S(6, 6);
	std::complex<double> S11 = Sss * (exp(1i * kVec[0] * aCC) + 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3.0) / 2.0 * kVec[1] * aCC));
	std::complex<double> S12 = Ssp * (-exp(1i * kVec[0] * aCC) + exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3.0) / 2.0 * kVec[1] * aCC));
	std::complex<double> S13 = Ssp * (-1i * sqrt(3.0) * exp(-1i * kVec[0] * aCC / 2.0) * sin(sqrt(3.0) / 2.0 * kVec[1] * aCC));
	std::complex<double> S21 = -S12;
	std::complex<double> S22 = -Ssigma * exp(1i * kVec[0] * aCC) + (3.0 * Spi - Ssigma) / 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3.0) / 2.0 * kVec[1] * aCC);
	std::complex<double> S23 = 1i * sqrt(3.0) / 2.0 * (Ssigma + Spi) * exp(-1i * kVec[0] * aCC / 2.0) * sin(sqrt(3.0) / 2.0 * kVec[1] * aCC);
	std::complex<double> S31 = -S13;
	std::complex<double> S32 = S23;
	std::complex<double> S33 = Spi * exp(1i * kVec[0] * aCC) + (Spi - 3.0 * Ssigma) / 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * kVec[1] * aCC);
	S << 1.0, 0.0, 0.0, S11, S12, S13,
		0.0, 1.0, 0.0, S21, S22, S23,
		0.0, 0.0, 1.0, S31, S32, S33,
		std::conj(S11), std::conj(S21), std::conj(S31), 1.0, 0.0, 0.0,
		std::conj(S12), std::conj(S22), std::conj(S32), 0.0, 1.0, 0.0,
		std::conj(S13), std::conj(S23), std::conj(S33), 0.0, 0.0, 1.0;
	return S;
}

Eigen::MatrixXcd Graphene::H_Matrix_Sigma(const Eigen::Vector3d& kVec) const
{
	Eigen::MatrixXcd H(6, 6);
	std::complex<double> H11 = Hss * (exp(1i * kVec[0] * aCC) + 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * kVec[1] * aCC));
	std::complex<double> H12 = Hsp * (-exp(1i * kVec[0] * aCC) + exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * kVec[1] * aCC));
	std::complex<double> H13 = Hsp * (-1i * sqrt(3.0) * exp(-1i * kVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * kVec[1] * aCC));
	std::complex<double> H21 = -H12;
	std::complex<double> H22 = -Hsigma * exp(1i * kVec[0] * aCC) + (3.0 * Hpi - Hsigma) / 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * kVec[1] * aCC);
	std::complex<double> H23 = 1i * sqrt(3.0) / 2.0 * (Hsigma + Hpi) * exp(-1i * kVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * kVec[1] * aCC);
	std::complex<double> H31 = -H13;
	std::complex<double> H32 = H23;
	std::complex<double> H33 = Hpi * exp(1i * kVec[0] * aCC) + (Hpi - 3.0 * Hsigma) / 2.0 * exp(-1i * kVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * kVec[1] * aCC);
	H << epsilon_2s, 0.0, 0.0, H11, H12, H13,
		0.0, epsilon_2p, 0.0, H21, H22, H23,
		0.0, 0.0, epsilon_2p, H31, H32, H33,
		std::conj(H11), std::conj(H21), std::conj(H31), epsilon_2s, 0.0, 0.0,
		std::conj(H12), std::conj(H22), std::conj(H32), 0.0, epsilon_2p, 0.0,
		std::conj(H13), std::conj(H23), std::conj(H33), 0.0, 0.0, epsilon_2p;
	return H;
}

Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> Graphene::EigenSolution_Pi(const Eigen::Vector3d& kVec, bool compute_eigenvectors) const
{
	Eigen::MatrixXcd H = H_Matrix_Pi(kVec);
	Eigen::MatrixXcd S = S_Matrix_Pi(kVec);
	return Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>(H, S, compute_eigenvectors);
}

Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> Graphene::EigenSolution_Sigma(const Eigen::Vector3d& kVec, bool compute_eigenvectors) const
{
	Eigen::MatrixXcd H = H_Matrix_Sigma(kVec);
	Eigen::MatrixXcd S = S_Matrix_Sigma(kVec);
	return Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>(H, S, compute_eigenvectors);
}

Graphene::Graphene(double workfunction)
{
	aCC			  = 1.42 * Angstrom;
	a			  = aCC * sqrt(3.0);
	b			  = 4.0 * M_PI / sqrt(3.0) / a;
	work_function = workfunction;
	N_cell		  = 0.5 * 5.0e25 / kg;

	lattice_vectors = {
		{3.0 / 2.0 * aCC, a / 2.0, 0.0},
		{3.0 / 2.0 * aCC, -a / 2.0, 0.0},
		{0.0, -a, 0.0}};
	reciprocal_lattice_vectors = {
		{2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / a, 0.0},
		{2.0 * M_PI / sqrt(3.0) / a, -2.0 * M_PI / a, 0.0}};
	nearest_neighbors = {
		{aCC, 0.0, 0.0},
		{-aCC / 2.0, a / 2.0, 0.0},
		{-aCC / 2.0, -a / 2.0, 0.0}};
	high_symmetry_point_G = {0.0, 0.0, 0.0};
	high_symmetry_point_M = {2.0 * M_PI / sqrt(3.0) / a, 0.0, 0.0};
	high_symmetry_point_K = {2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / 3.0 / a, 0.0};

	// Overlap and transfer integrals
	s		   = 0.129;
	sPrime	   = 0.00866032;
	Sss		   = 0.212;
	Ssp		   = 0.159931;
	Ssigma	   = 0.146;
	Spi		   = 0.129;
	t		   = -3.033 * eV;
	Hss		   = -6.769 * eV;
	Hsp		   = -5.580 * eV;
	Hsigma	   = -5.037 * eV;
	Hpi		   = -3.033 * eV;
	epsilon_2s = -8.868 * eV;
	epsilon_2p = 0.0;

	Zeff_2s		 = 4.59381;
	Zeff_2px_2py = 5.48626;
	Zeff_2pz	 = 4.02474;
}

std::vector<double> Graphene::Energy_Dispersion_Pi(const Eigen::Vector3d& kVec) const
{
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges = EigenSolution_Pi(kVec, false);
	std::vector<double> eigenvalues								   = {ges.eigenvalues()[0], ges.eigenvalues()[1]};
	return eigenvalues;
}

std::vector<double> Graphene::Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& kVec) const
{
	double fNorm = std::abs(f_aux(kVec));
	return {(epsilon_2p + t * fNorm) / (1.0 + s * fNorm), (epsilon_2p - t * fNorm) / (1.0 - s * fNorm)};
}

std::vector<double> Graphene::Energy_Dispersion_Sigma(const Eigen::Vector3d& kVec) const
{
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges = EigenSolution_Sigma(kVec, false);
	std::vector<double> eigenvalues								   = {ges.eigenvalues()[0], ges.eigenvalues()[1], ges.eigenvalues()[2], ges.eigenvalues()[3], ges.eigenvalues()[4], ges.eigenvalues()[5]};
	return eigenvalues;
}

double Graphene::Valence_Band_Energies(const Eigen::Vector3d& kVec, unsigned int energy_band)
{
	if(energy_band == 0)
		return Energy_Dispersion_Pi_Analytic(kVec)[0];
	else if(energy_band < 4)
		return Energy_Dispersion_Sigma(kVec)[energy_band - 1];
	else
	{
		std::cerr << "Error in Graphene::Valence_Band_Energies(): energy band " << energy_band << "is not a valence band." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

std::vector<std::vector<double>> Graphene::Energy_Bands(unsigned int k_points)
{
	double kMin			   = 0.0;
	double kMax			   = 2.0 * M_PI / 3.0 / a * (3.0 + sqrt(3.0));
	std::vector<double> ks = libphysica::Linear_Space(kMin, kMax, k_points);
	std::vector<std::vector<double>> energy_bands;
	for(auto& k : ks)
	{
		Eigen::Vector3d kVec		= Path_1BZ(k);
		std::vector<double> E_pi	= Energy_Dispersion_Pi_Analytic(kVec);
		std::vector<double> E_sigma = Energy_Dispersion_Sigma(kVec);
		std::vector<double> aux		= {k, E_pi[0], E_pi[1], E_sigma[0], E_sigma[1], E_sigma[2], E_sigma[3], E_sigma[4], E_sigma[5]};

		energy_bands.push_back(aux);
	}
	return energy_bands;
}

bool Graphene::In_1BZ(const Eigen::Vector3d& kVec)
{
	if(std::fabs(kVec[0]) > 4.0 * M_PI / sqrt(3.0) / a)
		return false;
	else if(std::fabs(kVec[1]) > 4.0 * M_PI / 3.0 / a - std::fabs(kVec[0]) / sqrt(3.0))
		return false;
	else
		return true;
}

Eigen::Vector3d Graphene::Find_1BZ_Vector(const Eigen::Vector3d& kVec)
{
	if(In_1BZ(kVec))
		return kVec;
	else
	{
		Eigen::Vector3d kVec_1BZ = kVec;
		kVec_1BZ[0]				 = std::fmod(kVec_1BZ[0], 2.0 * M_PI / sqrt(3.0) / a);
		kVec_1BZ[1]				 = std::fmod(kVec_1BZ[1], 4.0 * M_PI / 3.0 / a);
		if(In_1BZ(kVec_1BZ))
			return kVec_1BZ;
		else if(kVec_1BZ[0] < 0.0 && kVec_1BZ[1] < 0.0)
			kVec_1BZ += reciprocal_lattice_vectors[0];
		else if(kVec_1BZ[0] < 0.0 && kVec_1BZ[1] > 0.0)
			kVec_1BZ += reciprocal_lattice_vectors[1];
		else if(kVec_1BZ[0] > 0.0 && kVec_1BZ[1] < 0.0)
			kVec_1BZ -= reciprocal_lattice_vectors[1];
		else if(kVec_1BZ[0] > 0.0 && kVec_1BZ[1] > 0.0)
			kVec_1BZ -= reciprocal_lattice_vectors[0];
		return kVec_1BZ;
	}
}

std::complex<double> Graphene::Wavefunction_Momentum_Pi(const Eigen::Vector3d& lVec, const Eigen::Vector3d& kVec, int i) const
{
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges = EigenSolution_Pi(kVec);

	std::complex<double> C1 = ges.eigenvectors().col(i)[0];
	std::complex<double> C2 = ges.eigenvectors().col(i)[1];
	std::complex<double> f	= f_aux(kVec + lVec);

	double norm = 1.0;	 // GeneralizedSelfAdjointEigenSolver sets the eigenvectors such that C^* S C = 1
	return 1.0 / std::sqrt(norm) * (C1 + C2 * f) * Hydrogenic_Wavefunction_Momentum(lVec, "2pz", Zeff_2pz);
}

std::complex<double> Graphene::Wavefunction_Momentum_Pi_Analytic(const Eigen::Vector3d& lVec, const Eigen::Vector3d& kVec) const
{
	std::complex<double> f = f_aux(kVec);
	double phi			   = -atan2(f.imag(), f.real());

	f = f_aux(kVec + lVec);

	double norm = 1.0;
	for(int i = 0; i < 3; i++)
		norm += s * cos(phi + nearest_neighbors[i].dot(kVec));	 // + sPrime * cos(kVec.dot(lattice_vectors[i]));

	return pow(2.0 * norm, -0.5) * (1.0 + exp(1i * phi) * f) * Hydrogenic_Wavefunction_Momentum(lVec, "2pz", Zeff_2pz);
}

std::complex<double> Graphene::Wavefunction_Momentum_Sigma(const Eigen::Vector3d& lVec, const Eigen::Vector3d& kVec, int i) const
{

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges = EigenSolution_Sigma(kVec);

	std::complex<double> C1 = ges.eigenvectors().col(i)[0];
	std::complex<double> C2 = ges.eigenvectors().col(i)[1];
	std::complex<double> C3 = ges.eigenvectors().col(i)[2];
	std::complex<double> C4 = ges.eigenvectors().col(i)[3];
	std::complex<double> C5 = ges.eigenvectors().col(i)[4];
	std::complex<double> C6 = ges.eigenvectors().col(i)[5];
	std::complex<double> f	= f_aux(kVec + lVec);

	double norm = 1.0;	 // GeneralizedSelfAdjointEigenSolver sets the eigenvectors such that C^* S C = 1
	return 1.0 / std::sqrt(norm) * ((C1 + C4 * f) * Hydrogenic_Wavefunction_Momentum(lVec, "2s", Zeff_2s) + (C2 + C5 * f) * Hydrogenic_Wavefunction_Momentum(lVec, "2px", Zeff_2px_2py) + (C3 + C6 * f) * Hydrogenic_Wavefunction_Momentum(lVec, "2py", Zeff_2px_2py));
}

double Graphene::DM_Response(int band, const Eigen::Vector3d& qVec, const Eigen::Vector3d& kPrimeVec)
{
	std::complex<double> psi;

	Eigen::Vector3d lVec = qVec - kPrimeVec;

	// Determine the crystal momentum vector k
	Eigen::Vector3d kVec({kPrimeVec[0] - qVec[0], kPrimeVec[1] - qVec[1], 0.0});
	kVec = Find_1BZ_Vector(kVec);

	if(band == 0)
	{
		std::complex<double> f = f_aux(kVec);
		double phi_l		   = -atan2(f.imag(), f.real());
		double phi_2pz		   = Hydrogenic_Wavefunction_Momentum(lVec, "2pz", Zeff_2pz);
		double norm			   = 1.0;
		for(int i = 0; i < 3; i++)
			norm += s * cos(phi_l + nearest_neighbors[i].dot(kVec));
		double N_l = 1.0 / std::sqrt(norm);
		return std::pow(N_l * phi_2pz, 2.0) * (1.0 + cos(phi_l + nearest_neighbors[0].dot(kVec + lVec)));
	}
	else
	{
		int i														   = band - 1;
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges = EigenSolution_Sigma(kVec);

		std::complex<double> C1	 = ges.eigenvectors().col((i))[0];
		std::complex<double> C2	 = ges.eigenvectors().col((i))[1];
		std::complex<double> C3	 = ges.eigenvectors().col((i))[2];
		std::complex<double> C4	 = ges.eigenvectors().col((i))[3];
		std::complex<double> C5	 = ges.eigenvectors().col((i))[4];
		std::complex<double> C6	 = ges.eigenvectors().col((i))[5];
		Eigen::MatrixXcd S		 = S_Matrix_Sigma(kVec);
		double norm				 = ges.eigenvectors().col((i)).dot(S * ges.eigenvectors().col((i))).real();
		double N_l				 = 1.0 / std::sqrt(norm);
		std::complex<double> aux = std::exp(1i * nearest_neighbors[0].dot(kVec + lVec));
		psi						 = N_l * ((C1 + C4 * aux) * Hydrogenic_Wavefunction_Momentum(lVec, "2s", Zeff_2s) + (C2 + C5 * aux) * Hydrogenic_Wavefunction_Momentum(lVec, "2px", Zeff_2px_2py) + (C3 + C6 * aux) * Hydrogenic_Wavefunction_Momentum(lVec, "2py", Zeff_2px_2py));
		return std::norm(psi);
	}
}

}	// namespace graphene