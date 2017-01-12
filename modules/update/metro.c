
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
	sun_mat tmp1, tmp2, tmp3, tmp4, tmp5;
	sun_zero(*stap);
	for(unsigned int nu=0; nu<DIM; ++nu)
	{
		sun_zero(tmp1); sun_zero(tmp2); sun_zero(tmp3); sun_zero(tmp4); sun_zero(tmp5);
		if(nu != dir){
			//sun_mul_dag(tmp1, *pu[neib[n][dir]][nu], *pu[neib[n][nu]][dir]);
			//sun_mul_dag(tmp2, tmp1, *pu[n][nu]);
			//sun_self_add(*stap, tmp2);
			//sun_dag(tmp3, *pu[neib[neib[n][nu+DIM]][dir]][nu]);
			//sun_mul_dag(tmp4, tmp3, *pu[neib[n][nu+DIM]][dir]);
			//sun_mul(tmp5, tmp4, *pu[neib[n][nu+DIM]][nu]);
			//sun_self_add(*stap, tmp5);
			sun_mul(tmp1, *pu[n][nu], *pu[neib[n][nu]][dir]);
			sun_mul_dag(tmp2, tmp1, *pu[neib[n][dir]][nu]);
			sun_self_add(*stap, tmp2);
			sun_dag(tmp3, *pu[neib[n][nu+DIM]][nu]);
			sun_mul(tmp4, tmp3, *pu[neib[n][nu+DIM]][dir]);
			sun_mul(tmp5, tmp4, *pu[neib[neib[n][nu+DIM]][dir]][nu]);
			sun_self_add(*stap, tmp5);
		}
	}
}

static void epsball(sun_mat *u)
{
	sun_alg X;
	double *r;
	r = (double*)&X;
	double q[ALGVOL];
	ranlxd(q, ALGVOL);
	for(unsigned int i=0; i<ALGVOL; ++i)
	{
		r[i] = (2*q[i] - 1.0);
	}
	expx(runp.eps, &X, u);
}

double local_metr(int n, int dir, int m)
{
	int acc = 0;
	double p, R, trace, result;
	sun_mat deltaU, mult, tmp;
	sun_mat *stap;
	stap = malloc(sizeof(sun_mat));
	for(unsigned int i=0; i<m; ++i)
	{
		tmp = *pu[n][dir];
		epsball(&tmp);
		sun_sub(deltaU, tmp, *pu[n][dir]);
		staples(n, dir, stap);
		sun_mul_dag(mult, deltaU, *stap);
		sun_trace(trace, mult);
		R = exp((runp.beta/SUN)*trace);
		ranlxd(&p, 1);
		if(R >= 1 || p <= R)
		{
			*pu[n][dir] = tmp;
			++acc;
		}
	}
	free(stap);
	result = (double)acc/(double)m;
	return result;
}
#if(SUN==2)
void su2_relax(int n, int dir, su2mat *stap)
{
	double sqrtDet;
	su2mat stapDag, tmp, U, A, Adag;
	su2_dag(A, *stap);
	su2_det_sqrt(sqrtDet, A);
	su2_dble_mul(A, 1/sqrtDet);
	su2_dag(Adag, A);
	U = *pu[n][dir];
	su2_mat_mul_dag(tmp, Adag, U);
	su2_mat_mul(*stap, tmp, Adag);
}
#endif
/*void su3_relax(int n, int dir, su3mat *stap)
{

}*/
