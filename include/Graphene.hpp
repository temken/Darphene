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

	// Overlap and transfer integrals
	double s, Sss, Ssp, Ssigma, Spi, t, Hss, Hsp, Hsigma, Hpi, epsilon_2s, epsilon_2p;

	std::complex<double> Bloch_Wavefunction_A(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, double Zeff) const;
	std::complex<double> Bloch_Wavefunction_B(const Eigen::Vector3d& rVec, const Eigen::Vector3d& lVec, double Zeff) const;

	std::complex<double> f_aux(const Eigen::Vector3d& lVec) const;

	Eigen::MatrixXcd S_Matrix_Pi(const Eigen::Vector3d& lVec) const;
	Eigen::MatrixXcd H_Matrix_Pi(const Eigen::Vector3d& lVec) const;

  public:
	Graphene();

	std::vector<double> Energy_Dispersion_Pi(const Eigen::Vector3d& lVec) const;
	std::vector<double> Energy_Dispersion_Pi_Analytic(const Eigen::Vector3d& lVec) const;
};
}	// namespace graphene

#endif