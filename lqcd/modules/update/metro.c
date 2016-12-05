
/*******************************************************************************
*
* File metro.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to determine the exponential function exp(X),
* where X is an element of the Lie algebra of SU(N), correctly to O(X^2).
*
* The routines included are similar to those used by Martin Luescher in
* the DD-HMC code.
*
* Externally accessible functions:
*
*
*******************************************************************************/

#define METRO_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include"headers.h"
#include"modules.h"
#include "ranlxd.h"

void staples(int n, int dir, sun_mat *stap)
{
	sun_mat res, tmp1, tmp2, tmp3, tmp4, tmp5;
	sun_zero(res);
	for(unsigned int nu=0; nu<DIM; ++nu)
	{
		sun_zero(tmp1); sun_zero(tmp2); sun_zero(tmp3); sun_zero(tmp4); sun_zero(tmp5);
		if(nu == dir){
		}
		else{
			//+nu
			sun_mul_dag(tmp1, *pu[neib[n][dir]][nu], *pu[neib[n][nu]][dir]);
			sun_mul_dag(tmp2, tmp1, *pu[n][nu]);
			sun_self_add(res, tmp2);
			//-nu
			su3_dag(tmp3, *pu[neib[neib[n][nu+DIM]][dir]][nu]);
			sun_mul_dag(tmp4, tmp3, *pu[neib[n][nu+DIM]][dir]);
			sun_mul(tmp5, tmp4, *pu[neib[n][nu+DIM]][nu]);
			sun_self_add(res, tmp5);
		}
	}
	stap = &res;
}

static void epsball(sun_mat *u)
{
	//double eps = 0.1;
	sun_alg X;
	double *r;
	r = (double*)&X;
	double q[ALGVOL];
	ranlxd(q, ALGVOL);
	for(unsigned int i=0; i<ALGVOL; ++i)
	{
		r[i] = (2*q[i] - 1.0) * runp.eps;
	}
	expx(runp.eps, &X, u);
}

double local_metr(int n, int dir, int m)
{
	int acc = 0;
	double p, R, trace;
	sun_mat deltaU, stap, mult, tmp;
	for(unsigned int i=0; i<m; ++i)
	{
		tmp = *pu[n][dir];
		//logging("before = %f\n", tmp.c11.re);
		epsball(&tmp);
		//logging("after = %f\n", tmp.c11.re);
		sun_sub(deltaU, tmp, *pu[n][dir]);
		staples(n, dir, &stap);
		sun_mul_dag(mult, deltaU, stap);
		sun_trace(trace, mult);
		//logging("trace = %f\n", trace);
		R = exp((runp.beta/SUN)*trace);
		//logging("R = %f\n", R);
		ranlxd(&p, 1);
		if(R >= 1 || p <= R)
		{
			*pu[n][dir] = tmp;
			++acc;
		}
	}
	epsball(pu[0][0]);
	double result = (double)acc/(double)m;
	return result;
}
