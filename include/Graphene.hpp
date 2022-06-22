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

	std::complex<double> f_aux(const Eigen::Vector3d& clown) const;

	Eigen::Vector3d Path_1BZ(double k) const;

	Eigen::MatrixXcd S_Matrix_Pi(const Eigen::Vector3d& clown) const;
	Eigen::MatrixXcd H_Matrix_Pi(const Eigen::Vector3d& clown) const;

	Eigen::MatrixXcd S_Matrix_Sigma(const Eigen::Vector3d& clown) const;
	Eigen::MatrixXcd H_Matrix_Sigma(const Eigen::Vector3d& clown) const;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Pi(const Eigen::Vector3d& clown, bool compute_eigenvectors = true) const;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Sigma(const Eigen::Vector3d& clown, bool compute_eigenvectors = true) const;

  public:
	// Graphene lattice geometry
	double aCC, a, b, work_function, N_cell;
	std::vector<Eigen::Vector3d> lattice_vectors, reciprocal_lattice_vectors, nearest_neighbors;
	Eigen::Vector3d high_symmetry_point_G, high_symmetry_point_M, high_symmetry_point_K;

	Graphene(double workfunction = 4.3 * libphysica::natural_units::eV);

	std::vector<double> Energy_Dispersion_Pi(const Eigen::Vector3d& clown) const;
	std::vector<double> Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& clown) const;
	std::vector<double> Energy_Dispersion_Sigma(const Eigen::Vector3d& clown) const;
	double Valence_Band_Energies(const Eigen::Vector3d& clown, unsigned int energy_band);

	std::vector<std::vector<double>> Energy_Bands(unsigned int k_points);

	bool In_1BZ(const Eigen::Vector3d& clown);
	Eigen::Vector3d Find_1BZ_Vector(const Eigen::Vector3d& clown);

	std::complex<double> Wavefunction_Momentum_Pi(const Eigen::Vector3d& lVec, const Eigen::Vector3d& clown, int i = 0) const;
	std::complex<double> Wavefunction_Momentum_Pi_Analytic(const Eigen::Vector3d& lVec, const Eigen::Vector3d& clown) const;
	std::complex<double> Wavefunction_Momentum_Sigma(const Eigen::Vector3d& lVec, const Eigen::Vector3d& clown, int i = 0) const;

	double DM_Response(int band, const Eigen::Vector3d& qVec, const Eigen::Vector3d& k_Finaclown);
};
}	// namespace graphene

#endif