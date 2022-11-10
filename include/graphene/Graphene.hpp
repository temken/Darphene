#ifndef __Graphene_hpp_
#define __Graphene_hpp_

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <cmath>
#include <complex>

#include "libphysica/Natural_Units.hpp"

#include "Carbon_Wavefunctions.hpp"

namespace graphene
{

class Graphene
{
  protected:
	// Overlap and transfer integrals
	double s, sPrime, Sss, Ssp, Ssigma, Spi, t, Hss, Hsp, Hsigma, Hpi, epsilon_2s, epsilon_2p;

	std::complex<double> f_aux(const Eigen::Vector3d& kVec) const;

	Eigen::Vector3d Path_1BZ(double k) const;

	Eigen::MatrixXcd S_Matrix_Pi(const Eigen::Vector3d& kVec) const;
	Eigen::MatrixXcd H_Matrix_Pi(const Eigen::Vector3d& kVec) const;

	Eigen::MatrixXcd S_Matrix_Sigma(const Eigen::Vector3d& kVec) const;
	Eigen::MatrixXcd H_Matrix_Sigma(const Eigen::Vector3d& kVec) const;

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Pi(const Eigen::Vector3d& kVec, bool compute_eigenvectors = true) const;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> EigenSolution_Sigma(const Eigen::Vector3d& kVec, bool compute_eigenvectors = true) const;

	Carbon_Wavefunctions* carbon_wavefunctions = {nullptr};

	std::vector<double> normalization_corrections;
	void Compute_Normalization_Corrections();

  public:
	// Graphene lattice geometry
	double aCC, a, b, work_function, N_cell;
	std::vector<Eigen::Vector3d> lattice_vectors, reciprocal_lattice_vectors, nearest_neighbors;
	Eigen::Vector3d high_symmetry_point_G, high_symmetry_point_M, high_symmetry_point_K;

	Graphene(const std::string& wavefunctions, double workfunction = 4.3 * libphysica::natural_units::eV);

	double Overlap_Integral(const std::string& parameter);

	std::vector<double> Energy_Dispersion_Pi(const Eigen::Vector3d& kVec) const;
	std::vector<double> Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& kVec) const;
	std::vector<double> Energy_Dispersion_Sigma(const Eigen::Vector3d& kVec) const;
	double Valence_Band_Energies(const Eigen::Vector3d& kVec, unsigned int energy_band);

	std::vector<std::vector<double>> Energy_Bands(unsigned int k_points);

	bool In_1BZ(const Eigen::Vector3d& kVec);
	Eigen::Vector3d Find_G_Vector(const Eigen::Vector3d& kVec);

	double Material_Response_Function(int band, const Eigen::Vector3d& lVec);
	double Material_Response_Function(const Eigen::Vector3d& lVec);
};
}	// namespace graphene

#endif