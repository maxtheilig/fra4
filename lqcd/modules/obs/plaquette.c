
/*******************************************************************************
*
* File init.c
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
* double plaquette(void)
*      Computes the average plaquette.
*
* double gauge_action(void)
* 	   Computes the Wilson gauge action.
*
*******************************************************************************/

#define PLAQUETTE_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include"headers.h"
#include"modules.h"

double plaquette(void)
{
	//logging("computing the average plaquette\n");
	double result, factor;
	sun_mat tmp1, tmp2, tmp3, result_matrix;
	sun_zero(result_matrix);
	result = 0.0;
	for(unsigned int x = 0; x < VOL; x++) {
		for(unsigned int mu = 0; mu < DIM; mu++) {
			for(unsigned int nu = mu + 1; nu < DIM; nu++) {
				sun_mul(tmp1, *pu[x][mu], *pu[neib[x][mu]][nu]);
				sun_mul_dag(tmp2, tmp1, *pu[neib[x][nu]][mu]);
				sun_mul_dag(tmp3, tmp2, *pu[x][nu]);
				sun_self_add(result_matrix, tmp3);
			}
		}
	}
	sun_trace(result, result_matrix);
	factor = 1.0 / (NPLAQ * SUN * VOL);
	result *= factor;
	return result;
}

double gauge_action(void)
{
	double result;
	result = runp.beta * NPLAQ * VOL * ( 1 - plaquette());
	return result;
}
