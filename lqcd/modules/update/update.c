
/*******************************************************************************
*
* File update.c
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

#define UPDATE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include"headers.h"
#include"modules.h"
#include "ranlxd.h"

void update(int iup, int nup, int utype, int stype)
{
	double acc = 0.0;
	double plaq;
	for(unsigned int k=0; k<nup; ++k)
	{
		if(utype ==0)
		{
			for(unsigned int n=0; n<VOL; ++n)
			{
				for(unsigned int dir=0; dir<DIM; ++dir)
				{
					acc += local_metr(n, dir, iup);
				}
			}
		}
		plaq = plaquette();
		logging("average plaq = %f\n", plaq);
	}
	acc /= VOL*DIM*nup;
	logging("Acceptance Propability = %f\n", acc);
}
