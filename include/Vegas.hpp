#ifndef __Vegas_hpp_
#define __Vegas_hpp_

#include <functional>
#include <vector>

extern void vegas(std::vector<double>& regn, std::function<double(std::vector<double>&, const double)> fxn, const int init, const int ncall, const int itmx, const int nprn, double& tgral, double& sd, double& chi2a);
extern void brute_force(std::vector<double>& regn, std::function<double(std::vector<double>&, const double)> fxn, const int init, const int ncall, const int itmx, const int nprn, double& tgral, double& sd, double& chi2a);

#endif