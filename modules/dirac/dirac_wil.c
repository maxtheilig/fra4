
/*******************************************************************************
*
* File dirac_wil.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to initialise the fields and the error checks.
*
* Externally accessible functions:
*
*******************************************************************************/

#define DIRAC_WIL_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"headers.h"
#include"modules.h"

void apply_dirac_wil(sun_wferm *r, sun_wferm *s)
{
/*#if(DIM!=4)
	logging("wrong dimension for pseudofermions");
#else

#endif*/
	int dir;
	sun_wferm tmp1, tmp2, tmp3, tmp4, sNeib;
	sun_mat U, Uminus;
	double mEff;
	mEff = runp.m + 4;
	for(unsigned int x=0; x<VOL; ++x)
	{
		sunwferm_real_mult(r[x], mEff, s[x]);
		///////////////
		// mu = 0
		//////////////
		dir = 0;
		// mu = +0
		sNeib = s[neib[x][dir]];
		if(neib[x][dir]<x)
			sunwferm_real_mult(sNeib, -1.0, sNeib);
		U = *pu[x][dir];
		sunwferm_sun_mult(tmp1, U, sNeib);
		mul_sunwferm_mg0(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1-g0)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		// mu = -0
		sNeib = s[neib[x][dir+DIM]];
		if(neib[x][dir+DIM]>x)
			sunwferm_real_mult(sNeib, -1.0, sNeib);
		Uminus = *pu[neib[x][dir+DIM]][dir]; //U_(-dir)(x) = U_(dir)(x-dir)^dagger
		sunwferm_sun_dag_mult(tmp1, Uminus, sNeib);
		mul_sunwferm_g0(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1+g0)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		///////////////
		// mu = 1
		//////////////
		dir = 1;
		// mu = +1
		sNeib = s[neib[x][dir]];
		U = *pu[x][dir];
		sunwferm_sun_mult(tmp1, U, sNeib);
		mul_sunwferm_mg1(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1-g1)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		// mu = -1
		sNeib = s[neib[x][dir+DIM]];
		Uminus = *pu[neib[x][dir+DIM]][dir]; //U_(-dir)(x) = U_(dir)(x-dir)^dagger
		sunwferm_sun_dag_mult(tmp1, Uminus, sNeib);
		mul_sunwferm_g1(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1+g1)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		///////////////
		// mu = 2
		//////////////
		dir = 2;
		// mu = +2
		sNeib = s[neib[x][dir]];
		U = *pu[x][dir];
		sunwferm_sun_mult(tmp1, U, sNeib);
		mul_sunwferm_mg2(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1-g2)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		// mu = -2
		sNeib = s[neib[x][dir+DIM]];
		Uminus = *pu[neib[x][dir+DIM]][dir]; //U_(-dir)(x) = U_(dir)(x-dir)^dagger
		sunwferm_sun_dag_mult(tmp1, Uminus, sNeib);
		mul_sunwferm_g2(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1+g2)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		///////////////
		// mu = 3
		//////////////
		dir = 3;
		// mu = +3
		sNeib = s[neib[x][dir]];
		U = *pu[x][dir];
		sunwferm_sun_mult(tmp1, U, sNeib);
		mul_sunwferm_mg3(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1-g3)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
		// mu = -3
		sNeib = s[neib[x][dir+DIM]];
		Uminus = *pu[neib[x][dir+DIM]][dir]; //U_(-dir)(x) = U_(dir)(x-dir)^dagger
		sunwferm_sun_dag_mult(tmp1, Uminus, sNeib);
		mul_sunwferm_g3(tmp2, tmp1);
		sunwferm_add_single(tmp2, tmp1); //tmp2 = (1+g3)*U*sNeib
		sunwferm_real_mult(tmp3, 0.5, tmp2);
		sunwferm_sub_single(r[x], tmp3);
	}
}
