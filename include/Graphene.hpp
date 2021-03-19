#ifndef __Graphene_hpp_
#define __Graphene_hpp_

#include <Eigen/Geometry>
#include <cmath>
#include <complex>

// Headers from libphysica
#include "Natural_Units.hpp"

namespace graphene
{

class Graphene
{
  protected:
	// Graphene lattice geometry
	double aCC, a;
	std::vector<Eigen::Vector3d> lattice_vectors, reciprocal_lattice_vectors, nearest_neighbors;
	Eigen::Vector3d high_symmetry_point_G, high_symmetry_point_M, high_symmetry_point_K;

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

  public:
	Graphene();

	std::vector<double> Energy_Dispersion_Pi(const Eigen::Vector3d& lVec) const;
	std::vector<double> Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& lVec) const;

	std::vector<double> Energy_Dispersion_Sigma(const Eigen::Vector3d& lVec) const;

	std::vector<std::vector<double>> Energy_Bands(unsigned int k_points);

	std::complex<double> Wavefunction_Pi(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const;
	std::complex<double> Wavefunction_Pi_Analytic(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const;
	// std::complex<double> Wavefunction_Sigma(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec) const;
};
}	// namespace graphene

#endif