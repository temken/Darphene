#ifndef __Graphene_hpp_
#define __Graphene_hpp_

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <cmath>
#include <complex>

#include "libphysica/Natural_Units.hpp"

namespace graphene
{

class Graphene
{
  protected:
	// Overlap and transfer integrals
	double s, sPrime, Sss, Ssp, Ssigma, Spi, t, Hss, Hsp, Hsigma, Hpi, epsilon_2s, epsilon_2p;

	double Zeff_2s, Zeff_2px_2py, Zeff_2pz;
	std::complex<double> Bloch_Wavefunction_A(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, const std::string& orbital, double Zeff) const;
	std::complex<double> Bloch_Wavefunction_B(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, const std::string& orbital, double Zeff) const;

	std::complex<double> f_aux(const Eigen::Vector3d& lVec) const;

	Eigen::Vector3d Path_1BZ(double k) const;

	Eigen::MatrixXcd S_Matrix_Pi(const Eigen::Vector3d& lVec) const;
	Eigen::MatrixXcd H_Matrix_Pi(const Eigen::Vector3d& lVec) const;

	Eigen::MatrixXcd S_Matrix_Sigma(const Eigen::Vector3d& lVec) const;
	Eigen::MatrixXcd H_Matrix_Sigma(const Eigen::Vector3d& lVec) const;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Pi(const Eigen::Vector3d& lVec, bool compute_eigenvectors = true) const;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Sigma(const Eigen::Vector3d& lVec, bool compute_eigenvectors = true) const;

  public:
	// Graphene lattice geometry
	double aCC, a, b;
	std::vector<Eigen::Vector3d> lattice_vectors, reciprocal_lattice_vectors, nearest_neighbors;
	Eigen::Vector3d high_symmetry_point_G, high_symmetry_point_M, high_symmetry_point_K;

	Graphene();

	std::vector<double> Energy_Dispersion_Pi(const Eigen::Vector3d& lVec) const;
	std::vector<double> Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& lVec) const;
	std::vector<double> Energy_Dispersion_Sigma(const Eigen::Vector3d& lVec) const;
	double Valence_Band_Energies(const Eigen::Vector3d& lVec, unsigned int energy_band);

	std::vector<std::vector<double>> Energy_Bands(unsigned int k_points);

	bool In_1BZ(const Eigen::Vector3d& lVec);
	Eigen::Vector3d Find_1BZ_Vector(const Eigen::Vector3d& lVec);

	std::complex<double> Wavefunction_Pi(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, int i = 0) const;
	std::complex<double> Wavefunction_Pi_Analytic(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const;
	std::complex<double> Wavefunction_Sigma(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, int i = 0) const;

	std::complex<double> Wavefunction_Momentum_Pi(const Eigen::Vector3d& kVec, const Eigen::Vector3d& lVec, int i = 0) const;
	std::complex<double> Wavefunction_Momentum_Pi_Analytic(const Eigen::Vector3d& kVec, const Eigen::Vector3d& lVec) const;
	std::complex<double> Wavefunction_Momentum_Sigma(const Eigen::Vector3d& kVec, const Eigen::Vector3d& lVec, int i = 0) const;

	double DM_Response(int band, const Eigen::Vector3d& lVec, const Eigen::Vector3d& kVec);
	double DM_Response_corrected(int band, const Eigen::Vector3d& lVec, const Eigen::Vector3d& kVec);
};
}	// namespace graphene

#endif