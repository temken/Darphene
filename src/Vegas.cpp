#include "Vegas.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>

#include "Linear_Algebra.hpp"
#include "Statistics.hpp"

// Reference: These functions are taken from http://numerical.recipes/webnotes/nr3web9.pdf

// Utility routine used by vegas, to rebin a vector of densities contained in row j of xi into new bins defined by a vector r.
void rebin(const double rc, const int nd, std::vector<double>& r, std::vector<double>& xin, libphysica::Matrix& xi, const int j)
{
	int i, k = 0;
	double dr = 0.0, xn = 0.0, xo = 0.0;

	for(i = 0; i < nd - 1; i++)
	{
		while(rc > dr)
			dr += r[(++k) - 1];
		if(k > 1)
			xo = xi[j][k - 2];
		xn = xi[j][k - 1];
		dr -= rc;
		xin[i] = xn - (xn - xo) * dr / r[k - 1];
	}
	for(i = 0; i < nd - 1; i++)
		xi[j][i] = xin[i];
	xi[j][nd - 1] = 1.0;
}

void vegas(std::vector<double>& regn, std::function<double(std::vector<double>&, const double)> fxn, const int init, const int ncall, const int itmx, const int nprn, double& tgral, double& sd, double& chi2a)
{
	// Best make everything static, allowing restarts.
	static const int NDMX = 50, MXDIM = 10;
	static const double ALPH = 1.5, TINY = 1.0e-30;
	static int i, it, j, k, mds, nd, ndo, ng, npg;
	static double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti;
	static double tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt;
	static std::vector<int> ia(MXDIM), kg(MXDIM);
	static std::vector<double> dt(MXDIM), dx(MXDIM), r(NDMX), x(MXDIM), xin(NDMX);
	static libphysica::Matrix d(NDMX, MXDIM), di(NDMX, MXDIM), xi(MXDIM, NDMX);

	// Initialize  captive, static random number generator
	std::random_device rd;
	std::mt19937 PRNG(rd());

	int ndim = regn.size() / 2;
	if(init <= 0)
	{
		mds = ndo = 1;
		for(j = 0; j < ndim; j++)
			xi[j][0] = 1.0;
	}
	if(init <= 1)
		si = swgt = schi = 0.0;
	if(init <= 2)
	{
		nd = NDMX;
		ng = 1;
		if(mds != 0)
		{
			ng	= int(pow(ncall / 2.0 + 0.25, 1.0 / ndim));
			mds = 1;
			if((2 * ng - NDMX) >= 0)
			{
				mds = -1;
				npg = ng / NDMX + 1;
				nd	= ng / npg;
				ng	= npg * nd;
			}
		}
		for(k = 1, i = 0; i < ndim; i++)
			k *= ng;
		npg	  = std::max(int(ncall / k), 2);
		calls = double(npg) * double(k);
		dxg	  = 1.0 / ng;
		for(dv2g = 1, i = 0; i < ndim; i++)
			dv2g *= dxg;
		dv2g = calls * dv2g * calls * dv2g / npg / npg / (npg - 1.0);
		xnd	 = nd;
		dxg *= xnd;
		xjac = 1.0 / calls;
		for(j = 0; j < ndim; j++)
		{
			dx[j] = regn[j + ndim] - regn[j];
			xjac *= dx[j];
		}
		if(nd != ndo)
		{
			for(i = 0; i < std::max(nd, ndo); i++)
				r[i] = 1.0;
			for(j = 0; j < ndim; j++)
				rebin(ndo / xnd, nd, r, xin, xi, j);
			ndo = nd;
		}
		if(nprn >= 0)
		{
			std::cout << " Input parameters for vegas";
			std::cout << "  ndim= " << std::setw(4) << ndim;
			std::cout << "  ncall= " << std::setw(8) << calls << std::endl;
			std::cout << std::setw(34) << "  it=" << std::setw(5) << it;
			std::cout << "  itmx=" << std::setw(5) << itmx << std::endl;
			std::cout << std::setw(34) << "  nprn=" << std::setw(5) << nprn;
			std::cout << "  ALPH=" << std::setw(9) << ALPH << std::endl;
			std::cout << std::setw(34) << "  mds=" << std::setw(5) << mds;
			std::cout << "  nd=" << std::setw(5) << nd << std::endl;
			for(j = 0; j < ndim; j++)
			{
				std::cout << std::setw(30) << " x1[" << std::setw(2) << j;
				std::cout << "]= " << std::setw(11) << regn[j] << " xu[";
				std::cout << std::setw(2) << j << "]= ";
				std::cout << std::setw(11) << regn[j + ndim] << std::endl;
			}
		}
	}
	for(it = 0; it < itmx; it++)
	{
		ti = tsi = 0.0;
		for(j = 0; j < ndim; j++)
		{
			kg[j] = 1;
			for(i = 0; i < nd; i++)
				d[i][j] = di[i][j] = 0.0;
		}
		for(;;)
		{
			fb = f2b = 0.0;
			for(k = 0; k < npg; k++)
			{
				wgt = xjac;
				for(j = 0; j < ndim; j++)
				{
					xn	  = (kg[j] - libphysica::Sample_Uniform(PRNG)) * dxg + 1.0;
					ia[j] = std::max(std::min(int(xn), NDMX), 1);
					if(ia[j] > 1)
					{
						xo = xi[j][ia[j] - 1] - xi[j][ia[j] - 2];
						rc = xi[j][ia[j] - 2] + (xn - ia[j]) * xo;
					}
					else
					{
						xo = xi[j][ia[j] - 1];
						rc = (xn - ia[j]) * xo;
					}
					x[j] = regn[j] + rc * dx[j];
					wgt *= xo * xnd;
				}
				f  = wgt * fxn(x, wgt);
				f2 = f * f;
				fb += f;
				f2b += f2;
				for(j = 0; j < ndim; j++)
				{
					di[ia[j] - 1][j] += f;
					if(mds >= 0)
						d[ia[j] - 1][j] += f2;
				}
			}
			f2b = sqrt(f2b * npg);
			f2b = (f2b - fb) * (f2b + fb);
			if(f2b <= 0.0)
				f2b = TINY;
			ti += fb;
			tsi += f2b;
			if(mds < 0)
			{
				for(j = 0; j < ndim; j++)
					d[ia[j] - 1][j] += f2b;
			}
			for(k = ndim - 1; k >= 0; k--)
			{
				kg[k] %= ng;
				if(++kg[k] != 1)
					break;
			}
			if(k < 0)
				break;
		}
		tsi *= dv2g;
		wgt = 1.0 / tsi;
		si += wgt * ti;
		schi += wgt * ti * ti;
		swgt += wgt;
		tgral = si / swgt;
		chi2a = (schi - si * tgral) / (it + 0.0001);
		if(chi2a < 0.0)
			chi2a = 0.0;
		sd	= sqrt(1.0 / swgt);
		tsi = sqrt(tsi);
		if(nprn >= 0)
		{
			std::cout << " iteration no. " << std::setw(3) << (it + 1);
			std::cout << " : integral = " << std::setw(14) << ti;
			std::cout << " +/- " << std::setw(9) << tsi << std::endl;
			std::cout << " all iterations:  "
					  << " integral =";
			std::cout << std::setw(14) << tgral << "+-" << std::setw(9) << sd;
			std::cout << " chi**2/IT n =" << std::setw(9) << chi2a << std::endl;
			if(nprn != 0)
			{
				for(j = 0; j < ndim; j++)
				{
					std::cout << " DATA FOR axis  " << std::setw(2) << j << std::endl;
					std::cout << "     X      delta i          X      delta i";
					std::cout << "          X       deltai" << std::endl;
					for(i = nprn / 2; i < nd - 2; i += nprn + 2)
					{
						std::cout << std::setw(8) << xi[j][i] << std::setw(12) << di[i][j];
						std::cout << std::setw(12) << xi[j][i + 1] << std::setw(12) << di[i + 1][j];
						std::cout << std::setw(12) << xi[j][i + 2] << std::setw(12) << di[i + 2][j];
						std::cout << std::endl;
					}
				}
			}
		}
		for(j = 0; j < ndim; j++)
		{
			xo		= d[0][j];
			xn		= d[1][j];
			d[0][j] = (xo + xn) / 2.0;
			dt[j]	= d[0][j];
			for(i = 2; i < nd; i++)
			{
				rc			= xo + xn;
				xo			= xn;
				xn			= d[i][j];
				d[i - 1][j] = (rc + xn) / 3.0;
				dt[j] += d[i - 1][j];
			}
			d[nd - 1][j] = (xo + xn) / 2.0;
			dt[j] += d[nd - 1][j];
		}
		for(j = 0; j < ndim; j++)
		{
			rc = 0.0;
			for(i = 0; i < nd; i++)
			{
				if(d[i][j] < TINY)
					d[i][j] = TINY;
				r[i] = pow((1.0 - d[i][j] / dt[j]) /
							   (log(dt[j]) - log(d[i][j])),
						   ALPH);
				rc += r[i];
			}
			rebin(rc / xnd, nd, r, xin, xi, j);
		}
	}
}

void brute_force(std::vector<double>& regn, std::function<double(std::vector<double>&, const double)> fxn, const int init, const int ncall, const int itmx, const int nprn, double& tgral, double& sd, double& chi2a)
{
	int dim = regn.size() / 2.0;
	std::random_device rd;
	std::mt19937 PRNG(rd());

	double volume = 1.0;
	for(int i = 0; i < dim; i++)
		volume *= (regn[i + dim] - regn[i]);
	std::cout << "Volume (bf) = " << volume << std::endl;

	double sum	 = 0.0;
	double sum_2 = 0.0;
	for(int i = 0; i < ncall; i++)
	{
		std::vector<double> args(dim);
		for(int j = 0; j < dim; j++)
			args[j] = regn[j] + libphysica::Sample_Uniform(PRNG) * (regn[j + dim] - regn[j]);
		double fct = fxn(args, 0.0);
		sum += volume * fct;
		sum_2 += volume * volume * fct * fct;
	}
	tgral = sum / ncall;
	sd	  = sqrt((sum_2 / ncall - tgral * tgral) / (ncall - 1.0));
	chi2a = 0.0;
}