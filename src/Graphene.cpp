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
	return Hydrogenic_Wavefunction_2s(rVec, Zeff);
}

std::complex<double> Graphene::Bloch_Wavefunction_B(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, double Zeff) const
{
	std::complex<double> Phi = 0.0;
	for(auto& R : nearest_neighbors)
		Phi += exp(1i * lVec.dot(R)) * Hydrogenic_Wavefunction_2s(rVec - R, Zeff);
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

Graphene::Graphene()
{
	aCC = 1.42 * Angstrom;
	a	= aCC * sqrt(3.0);

	lattice_vectors = {
		{3.0 / 2.0 * aCC, a / 2.0, 0.0},
		{3.0 / 2.0 * aCC, -a / 2.0, 0.0}};
	reciprocal_lattice_vectors = {
		{2.0 * M_PI / sqrt(3.0) / a, 2.0 * M_PI / a, 0.0},
		{2.0 * M_PI / sqrt(3.0) / a, -2.0 * M_PI / a, 0.0}};
	nearest_neighbors = {
		{aCC, 0.0, 0.0},
		{-aCC / 2.0, a / 2.0, 0.0},
		{-aCC / 2.0, -a / 2.0, 0.0}};

	// Overlap and transfer integrals
	s		   = 0.129;
	Sss		   = 0.212;
	Ssp		   = 0.102;
	Ssigma	   = 0.146;
	Spi		   = 0.129;
	t		   = -3.033 * eV;
	Hss		   = -6.769 * eV;
	Hsp		   = -5.580 * eV;
	Hsigma	   = -5.037 * eV;
	Hpi		   = -3.033 * eV;
	epsilon_2s = -8.868 * eV;
	epsilon_2p = 0;
}

std::vector<double> Graphene::Energy_Dispersion_Pi(const Eigen::Vector3d& lVec) const
{
	Eigen::MatrixXcd S = S_Matrix_Pi(lVec);
	Eigen::MatrixXcd H = H_Matrix_Pi(lVec);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(H, S, false);
	std::vector<double> eigenvalues = {ges.eigenvalues()[0], ges.eigenvalues()[1]};
	std::sort(eigenvalues.begin(), eigenvalues.end());

	return eigenvalues;
}

std::vector<double> Graphene::Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& lVec) const
{
	double fNorm = std::abs(f_aux(lVec));
	return {(epsilon_2p + t * fNorm) / (1.0 + s * fNorm), (epsilon_2p - t * fNorm) / (1.0 - s * fNorm)};
}

}	// namespace graphene