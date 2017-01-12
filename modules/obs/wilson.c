
/*******************************************************************************
*
* File wilson.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to measure observables.
*
* Externally accessible functions:
*
*******************************************************************************/

#define WILSON_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include"headers.h"
#include"modules.h"

static sun_mat *iu1[VOL][DIM], *iu2[VOL][DIM];

static double compute_wloop(int n, int d1, int d2, int l1, int l2)
{
	double res;
	sun_mat S1, S2, T1, T2, tmp, tmp1;
	sun_unit(S1); sun_unit(S2); sun_unit(T1); sun_unit(T2);
	for(unsigned int k=0; k<l2; ++k)
	{
		sun_mul(tmp, S1, *iu2[n][d2]);
		S1 = tmp;
		n = neib[n][d2];
	}
	for(unsigned int k=0; k<l1; ++k)
	{
		sun_mul(tmp, T1, *iu1[n][d1]);
		T1 = tmp;
		n = neib[n][d1];
	}
	for(unsigned int k=0; k<l2; ++k)
	{
		n = neib[n][d2+DIM];
		sun_mul_dag(tmp, S2, *iu2[n][d2]);
		S2 = tmp;
	}
	n = neib[n][d1+DIM];
	for(unsigned int k=0; k<l1; ++k)
	{
		n = neib[n][d1+DIM];
		sun_mul_dag(tmp, T2, *iu1[n][d1]);
		T2 = tmp;
	}
	sun_mul(tmp, S1, T1);
	sun_mul(tmp1, tmp, S2);
	sun_mul(tmp, tmp1, T2);
	sun_trace(res, tmp);
	res /= SUN;
	return res;
}

void meas_wilson(double *w)
{
	int tDir = 0; int sDir = 3;
	int n;
	int i = 0;
	allocate_gauge(iu1); allocate_gauge(iu2);
	copy_gauge_field(pu, iu1); copy_gauge_field(pu, iu2);
	//smearing_APE_temporal(inlmeas.ismt, inlmeas.smpart, iu1);
	//smearing_APE_spatial(inlmeas.isms, inlmeas.smpars, iu2);
	for(unsigned int tempExt=inlmeas.ts; tempExt<=inlmeas.tf; tempExt+=inlmeas.dt){
	for(unsigned int spatExt=inlmeas.rs; spatExt<=inlmeas.rf; spatExt+=inlmeas.dr){
		//logging("temporalExtend = %i \t spatialExtend = %i\n", tempExt, spatExt);
		w[i] = 0.0;
		for(unsigned int t=0; t<(LENGT-tempExt+1); t+=tempExt)
		{
			for(unsigned int x=0; x<(lat.barea[0]-4); ++x)
			{
				n = x + t * lat.barea[0];
				w[i] += compute_wloop(n, tDir, sDir, tempExt, spatExt);
			}
		}
		w[i] /= lat.barea[0] * (LENGT-tempExt)/tempExt;
		logging("W = %f\n", w[i]);
		++i;
	}
	}

	free_gauge(iu1); free_gauge(iu2);
}
