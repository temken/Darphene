#ifndef __Carbon_Wavefunctions_hpp_
#define __Carbon_Wavefunctions_hpp_

#include <Eigen/Geometry>
#include <complex>
#include <string>

namespace graphene
{

// 1. Base class for atomic wavefunctions of carbon
class Carbon_Wavefunctions
{
  protected:
  public:
	Carbon_Wavefunctions() {};

	// Position space wavefunctions
	virtual double Wavefunction_Position_2s(const Eigen::Vector3d& xVec) { return 0; };
	virtual double Wavefunction_Position_2px(const Eigen::Vector3d& xVec) { return 0; };
	virtual double Wavefunction_Position_2py(const Eigen::Vector3d& xVec) { return 0; };
	virtual double Wavefunction_Position_2pz(const Eigen::Vector3d& xVec) { return 0; };

	// Momentum space wavefunctions
	virtual std::complex<double> Wavefunction_Momentum_2s(const Eigen::Vector3d& kVec) { return 0; };
	virtual std::complex<double> Wavefunction_Momentum_2px(const Eigen::Vector3d& kVec) { return 0; };
	virtual std::complex<double> Wavefunction_Momentum_2py(const Eigen::Vector3d& kVec) { return 0; };
	virtual std::complex<double> Wavefunction_Momentum_2pz(const Eigen::Vector3d& kVec) { return 0; };
};

// 2. Hydrogenic wavefunctions
class Hydrogenic : public Carbon_Wavefunctions
{
  protected:
	double Zeff_2s, Zeff_2pxpy, Zeff_2pz;

  public:
	Hydrogenic(double Z2s, double Z2pxpy, double Z2pz)
	: Zeff_2s(Z2s), Zeff_2pxpy(Z2pxpy), Zeff_2pz(Z2pz) {};

	// Position space wavefunctions
	virtual double Wavefunction_Position_2s(const Eigen::Vector3d& xVec) override;
	virtual double Wavefunction_Position_2px(const Eigen::Vector3d& xVec) override;
	virtual double Wavefunction_Position_2py(const Eigen::Vector3d& xVec) override;
	virtual double Wavefunction_Position_2pz(const Eigen::Vector3d& xVec) override;

	// Momentum space wavefunctions
	virtual std::complex<double> Wavefunction_Momentum_2s(const Eigen::Vector3d& kVec) override;
	virtual std::complex<double> Wavefunction_Momentum_2px(const Eigen::Vector3d& kVec) override;
	virtual std::complex<double> Wavefunction_Momentum_2py(const Eigen::Vector3d& kVec) override;
	virtual std::complex<double> Wavefunction_Momentum_2pz(const Eigen::Vector3d& kVec) override;
};

}	// namespace graphene

#endif