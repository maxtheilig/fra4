
/*******************************************************************************
*
* File update.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Externally accessible functions:
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
	double r[2];
	int x,y;
	for(unsigned int k=0; k<nup; ++k)
	{
		if(utype == 0)
		{
			for(unsigned int n=0; n<VOL; ++n)
			{
				for(unsigned int dir=0; dir<DIM; ++dir)
				{
					acc += local_metr(n, dir, 1);
				}
			}
		}
		if(utype == 1)
		{
			for(unsigned int i=0; i<(VOL*DIM); ++i)
			{
				ranlxd(r, 2);
				x = (int) r[0] * VOL;
				y = (int) r[1] * DIM;
				acc += local_metr(x, y ,1);
			}
		}
	}
	acc /= VOL*DIM*nup;
	//logging("Acceptance Propability = %f\n", acc);
}
