#include "Graphene.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>

#include "Hydrogenic_Wavefunctions.hpp"

namespace graphene
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

std::complex<double> Graphene::Bloch_Wavefunction_A(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, double Zeff) const
{
	return Hydrogenic_Wavefunction_2pz(rVec, Zeff);
}

std::complex<double> Graphene::Bloch_Wavefunction_B(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, double Zeff) const
{
	std::complex<double> Phi = 0.0;
	for(auto& R : nearest_neighbors)
		Phi += exp(1i * lVec.dot(R)) * Hydrogenic_Wavefunction_2pz(rVec - R, Zeff);
	return Phi;
}

std::complex<double> Graphene::f_aux(const Eigen::Vector3d& lVec) const
{
	return exp(1i * lVec[0] * a / sqrt(3.0)) + 2.0 * exp(-1i * lVec[0] * a / 2.0 / sqrt(3.0)) * cos(a * lVec[1] / 2.0);
}

Eigen::MatrixXcd Graphene::S_Matrix_Pi(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S(2, 2);
	S << 1.0, s * f_aux(lVec), s * std::conj(f_aux(lVec)), 1.0;
	return S;
}

Eigen::MatrixXcd Graphene::H_Matrix_Pi(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd H(2, 2);
	H << epsilon_2p, t * f_aux(lVec), t * std::conj(f_aux(lVec)), epsilon_2p;
	return H;
}

Eigen::MatrixXcd Graphene::S_Matrix_Sigma(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S(6, 6);
	std::complex<double> S11 = Sss * (exp(1i * lVec[0] * aCC) + 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> S12 = Ssp * (-exp(1i * lVec[0] * aCC) + exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> S13 = Ssp * (-1i * sqrt(3.0) * exp(-1i * lVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> S21 = -S12;
	std::complex<double> S22 = -Ssigma * exp(1i * lVec[0] * aCC) + (3.0 * Spi - Ssigma) / 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC);
	std::complex<double> S23 = 1i * sqrt(3.0) / 2.0 * (Ssigma + Spi) * exp(-1i * lVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * lVec[1] * aCC);
	std::complex<double> S31 = -S13;
	std::complex<double> S32 = S23;
	std::complex<double> S33 = Spi * exp(1i * lVec[0] * aCC) + (Spi - 3.0 * Ssigma) / 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC);
	S << 1.0, 0.0, 0.0, S11, S12, S13,
		0.0, 1.0, 0.0, S21, S22, S23,
		0.0, 0.0, 1.0, S31, S32, S33,
		std::conj(S11), std::conj(S21), std::conj(S31), 1.0, 0.0, 0.0,
		std::conj(S12), std::conj(S22), std::conj(S32), 0.0, 1.0, 0.0,
		std::conj(S13), std::conj(S23), std::conj(S33), 0.0, 0.0, 1.0;
	return S;
}

Eigen::MatrixXcd Graphene::H_Matrix_Sigma(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd H(6, 6);
	std::complex<double> H11 = Hss * (exp(1i * lVec[0] * aCC) + 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> H12 = Hsp * (-exp(1i * lVec[0] * aCC) + exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> H13 = Hsp * (-1i * sqrt(3.0) * exp(-1i * lVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * lVec[1] * aCC));
	std::complex<double> H21 = -H12;
	std::complex<double> H22 = -Hsigma * exp(1i * lVec[0] * aCC) + (3.0 * Hpi - Hsigma) / 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC);
	std::complex<double> H23 = 1i * sqrt(3.0) / 2.0 * (Hsigma + Hpi) * exp(-1i * lVec[0] * aCC / 2.0) * sin(sqrt(3) / 2.0 * lVec[1] * aCC);
	std::complex<double> H31 = -H13;
	std::complex<double> H32 = H23;
	std::complex<double> H33 = Hpi * exp(1i * lVec[0] * aCC) + (Hpi - 3.0 * Hsigma) / 2.0 * exp(-1i * lVec[0] * aCC / 2.0) * cos(sqrt(3) / 2.0 * lVec[1] * aCC);
	H << epsilon_2s, 0.0, 0.0, H11, H12, H13,
		0.0, epsilon_2p, 0.0, H21, H22, H23,
		0.0, 0.0, epsilon_2p, H31, H32, H33,
		std::conj(H11), std::conj(H21), std::conj(H31), epsilon_2s, 0.0, 0.0,
		std::conj(H12), std::conj(H22), std::conj(H32), 0.0, epsilon_2p, 0.0,
		std::conj(H13), std::conj(H23), std::conj(H33), 0.0, 0.0, epsilon_2p;
	return H;
}

Graphene::Graphene()
{
	aCC = 1.42 * Angstrom;
	a	= aCC * sqrt(3.0);

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
	epsilon_2p = 0;

	Zeff_2s		 = 4.59381;
	Zeff_2px_2py = 5.48626;
	Zeff_2pz	 = 4.02474;
}

std::vector<double> Graphene::Energy_Dispersion_Pi(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S = S_Matrix_Pi(lVec);
	Eigen::MatrixXcd H = H_Matrix_Pi(lVec);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(H, S, false);
	std::vector<double> eigenvalues = {ges.eigenvalues()[0], ges.eigenvalues()[1]};
	// std::sort(eigenvalues.begin(), eigenvalues.end());

	return eigenvalues;
}

std::vector<double> Graphene::Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& lVec) const
{
	double fNorm = std::abs(f_aux(lVec));
	return {(epsilon_2p + t * fNorm) / (1.0 + s * fNorm), (epsilon_2p - t * fNorm) / (1.0 - s * fNorm)};
}

std::vector<double> Graphene::Energy_Dispersion_Sigma(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S = S_Matrix_Sigma(lVec);
	Eigen::MatrixXcd H = H_Matrix_Sigma(lVec);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(H, S, false);
	std::vector<double> eigenvalues = {ges.eigenvalues()[0], ges.eigenvalues()[1], ges.eigenvalues()[2], ges.eigenvalues()[3], ges.eigenvalues()[4], ges.eigenvalues()[5]};
	// std::sort(eigenvalues.begin(), eigenvalues.end());

	return eigenvalues;
}

std::complex<double> Graphene::Wavefunction_Pi(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S = S_Matrix_Pi(lVec);
	Eigen::MatrixXcd H = H_Matrix_Pi(lVec);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(H, S);
	std::complex<double> C1 = ges.eigenvectors().col(0)[0];
	std::complex<double> C2 = ges.eigenvectors().col(0)[1];
	double normalization	= 1.0;
	return normalization * (C1 * Bloch_Wavefunction_A(rVec, lVec, Zeff_2pz) + C2 * Bloch_Wavefunction_B(rVec, lVec, Zeff_2pz));
}

std::complex<double> Graphene::Wavefunction_Pi_Analytic(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const
{
	std::complex<double> f = f_aux(lVec);
	double phi			   = -atan2(f.imag(), f.real());
	// double norm			   = 1.0;
	double norm = 2.0;
	for(int i = 0; i < 3; i++)
		norm += s * cos(phi + nearest_neighbors[i].dot(lVec)) + sPrime * cos(lVec.dot(lattice_vectors[i]));
	return pow(2.0 * norm, -0.5) * (Bloch_Wavefunction_A(rVec, lVec, Zeff_2pz) + exp(1i * phi) * Bloch_Wavefunction_B(rVec, lVec, Zeff_2pz));
}

}	// namespace graphene