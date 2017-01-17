
/*******************************************************************************
*
* File solv_cg.c
*
* Copyright (C) 2017 Bastian Brandt, Tim Breitenfelder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the CG solver for inversion of fermion dirac operators.
*
* Externally accessible functions:
*
* cg(...)
*
*
*******************************************************************************/
#define SOLV_CG_C

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"ranlxd.h"
#include"headers.h"
#include"modules.h"
#include"float.h"

int cg(sun_wferm *x,
       void (*A)(sun_wferm *r,sun_wferm *s),
	   void (*Ad)(sun_wferm *r,sun_wferm *s),sun_wferm *b,
	   double eps, int nmax)
{
	sun_wferm *rho, *p, *z, *tmp;
	unsigned int counter = 0;
	double delta, alpha, beta, zeta, zeta_tmp;

	alloc_wferm(&rho); alloc_wferm(&p); alloc_wferm(&z); alloc_wferm(&tmp);

	// b = D^dag * eta

	A(tmp, x); Ad(tmp, tmp);
	sunwferm_sub(*rho, *b, *tmp);				// rho_0 = eta - A * phi_0
	p = rho;									// p_0 = rho_0
	zeta = square_norm(rho);					// zeta = (rho_0)^dag * rho_0
	delta = eps * square_norm(b);				// delta = epsilon * ||eta||

	while(counter < nmax)
	{
		A(tmp, p); Ad(z, tmp);					// z = A * p_(k-1)
		alpha = zeta / scalar_prod(p, z).re;	// alpha_(k-1) = zeta / [(p_(k-1))^dag * z]
		sunwferm_real_mult(*tmp, alpha, *p);
		sunwferm_add(*x, *x, *tmp);				// phi_k = phi_(k-1) + alpha_(k-1) * p_(k-1)
		sunwferm_real_mult(*tmp, alpha, *z);
		sunwferm_sub(*rho, *rho, *tmp);			// rho_k = rho_(k-1) - alpha_(k-1) * z

		if(square_norm(rho) < delta) return 0;	// if |rho_k| < delta -> leave

		if(100 * square_norm(x) * DBL_EPSILON > delta) return -2; // inversion unfeasible

		zeta_tmp = square_norm(rho);			// zeta' = (rho_k)^dag * rho_k
		beta = zeta_tmp / zeta;					// beta_(k-1) = zeta' / zeta
		sunwferm_real_mult(*tmp, beta, *p);
		sunwferm_add(*p, *rho, *tmp);				// p_k = rho_k + beta_(k-1) * p_(k-1);
		zeta = zeta_tmp;						// zeta = zeta'

		counter++;
	}

	return -1;									// nmax reached

	free_wferm(&rho); free_wferm(&p); free_wferm(&z); free_wferm(&tmp);
}
