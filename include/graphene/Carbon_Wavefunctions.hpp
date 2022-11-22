#ifndef __Carbon_Wavefunctions_hpp_
#define __Carbon_Wavefunctions_hpp_

#include <Eigen/Geometry>
#include <complex>
#include <string>
#include <vector>

namespace graphene
{

// 1. Base class for atomic wavefunctions of carbon
class Carbon_Wavefunctions
{
  protected:
  public:
	std::string name;
	Carbon_Wavefunctions() {};

	// Position space wavefunctions
	virtual double Wavefunction_Position(const Eigen::Vector3d& xVec, const std::string& orbital) { return 0; };

	// Momentum space wavefunctions
	virtual std::complex<double> Wavefunction_Momentum(const Eigen::Vector3d& kVec, const std::string& orbital) { return 0; };

	// Normalization
	double Normalization_Position(const std::string& orbital);
	double Normalization_Momentum(const std::string& orbital);
};

// 2. Hydrogenic wavefunctions
class Hydrogenic : public Carbon_Wavefunctions
{
  protected:
	double Zeff_2s, Zeff_2pxpy, Zeff_2pz;

  public:
	Hydrogenic(double Z2s, double Z2pxpy, double Z2pz)
	: Zeff_2s(Z2s), Zeff_2pxpy(Z2pxpy), Zeff_2pz(Z2pz) { name = "Hydrogenic"; };

	// Position space wavefunctions
	virtual double Wavefunction_Position(const Eigen::Vector3d& xVec, const std::string& orbital) override;

	// Momentum space wavefunctions
	virtual std::complex<double> Wavefunction_Momentum(const Eigen::Vector3d& kVec, const std::string& orbital) override;
};

// 3. Roothaan-Hartree-Fock wavefunctions
extern std::complex<double> Hypergeometric_2F1(double a, double b, double c, double z);
extern double Slate_Type_Orbital_Position(double r, double Z, double n);
extern std::complex<double> Slate_Type_Orbital_Momentum(double k, double Z, double n, int l);

class Roothaan_Hartree_Fock : public Carbon_Wavefunctions
{
  protected:
	std::vector<std::vector<double>> C_nlj, Z_lj, n_lj;

  public:
	void Import_RHF_Coefficients();

	double Radial_Wavefunction_Position(double r, const std::string& orbital) const;
	std::complex<double> Radial_Wavefunction_Momentum(double p, const std::string& orbital) const;

  public:
	Roothaan_Hartree_Fock();

	// Position space wavefunctions
	virtual double Wavefunction_Position(const Eigen::Vector3d& xVec, const std::string& orbital) override;

	// Momentum space wavefunctions
	virtual std::complex<double> Wavefunction_Momentum(const Eigen::Vector3d& kVec, const std::string& orbital) override;
};

}	// namespace graphene

#endif