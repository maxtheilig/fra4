
/*******************************************************************************
*
* File heat.c
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

#define HEAT_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include"headers.h"
#include"modules.h"
#include"ranlxd.h"

void su2_update_heat(su2mat *u0)
{
	double x0_new, betaEff, r1[4], r2[2], alpha, sqrtDet, sqrtOneMinusX0Squared;
	su2mat stapleDagger, X, A, ADagger;
	su2_dag(stapleDagger, *u0);
	su2_det_sqrt(sqrtDet, stapleDagger);
	//logging("sqrtDet = %f\n", sqrtDet);
	betaEff = 2.0 * runp.beta / SUN;
	alpha = betaEff * sqrtDet;
	//logging("alpha = %f, betaEff = %f\n", alpha, betaEff);

	A = stapleDagger;
	su2_dble_mul(A, 1/sqrtDet);
	su2_dag(ADagger, A);

	int acc = 0;
	while(!acc)
	{
		//logging("chilling in while loop\n");
		ranlxd(r1, 4);
		x0_new = 1 + 1/alpha * ( log(1-r1[0]) + log(1-r1[1]) * pow(cos(2*PI*r1[2]),2) );
		//logging("x0_new = %f\n", x0_new);
		if(2*pow(r1[3],2) <= 1+x0_new) {X.c0 = x0_new; ++acc;}
	}

	ranlxd(r2, 2);
	sqrtOneMinusX0Squared = sqrt(1-pow(X.c0,2));
	X.c1 = sqrtOneMinusX0Squared * (1-2*r2[0]);
	X.c2 = sqrtOneMinusX0Squared * sqrt(r2[0]) * cos(2*PI*r2[1]);
	X.c3 = sqrtOneMinusX0Squared * sqrt(r2[0]) * sin(2*PI*r2[1]);

	su2_mat_mul(*u0, ADagger, X);
}

#if(SUN==3)
void su3_update_heat(int n, int dir, su3mat *stap)
{
	double sqrtDet;
	su3mat Z, U, U1, U2, Z1, Z2, W12_3x3, W23_3x3, W13_3x3;
	su2mat W12, W23, W13, Z12dag, Z23dag, Z13dag, A, *X, Z12, Z23, Z13;
	X = malloc(sizeof(su2mat));

	//get link variable
	U = *pu[n][dir];

	//compute Z = U * staple^+
	su3_mat_mul_dag(Z, U, *stap);

	//////////////////////////
	// subgroup (12)
	//////////////////////////
	extr_sub1_dag(Z12dag, Z);

	//compute A
	su2_dag(Z12, Z12dag);
	su2_det_sqrt(sqrtDet, Z12);
	A = Z12;
	su2_dble_mul(A, 1/sqrtDet);

	//compute X via su2 heatbath
	extr_sub1(*X, *stap);
	su2_update_heat(X);

	//compute W12(2x2)
	su2_mat_mul_dag(W12, *X, A);

	//compute U1 = W12 * U1, Z1 = W12 * Z
	create_su3_sub1(W12_3x3, W12);
	su3_mat_mul(U1, W12_3x3, U);
	su3_mat_mul(Z1, W12_3x3, Z);

	//////////////////////////
	// subgroup (23)
	//////////////////////////
	extr_sub2_dag(Z23dag, Z1);

	//compute A
	su2_dag(Z23, Z23dag);
	su2_det_sqrt(sqrtDet, Z23);
	A = Z23;
	su2_dble_mul(A, 1/sqrtDet);

	//compute X via su2 heatbath
	extr_sub2(*X, *stap);
	su2_update_heat(X);

	//compute W23(2x2)
	su2_mat_mul_dag(W23, *X, A);

	//compute U2 = W23 * U1, Z2 = W23 * Z1
	create_su3_sub2(W23_3x3, W23);
	su3_mat_mul(U2, W23_3x3, U1);
	su3_mat_mul(Z2, W23_3x3, Z1);

	//////////////////////////
	// subgroup (13)
	//////////////////////////
	extr_sub3_dag(Z13dag, Z2);

	//compute A
	su2_dag(Z13, Z13dag);
	su2_det_sqrt(sqrtDet, Z13);
	A = Z13;
	su2_dble_mul(A, 1/sqrtDet);

	//compute X via su2 heatbath
	extr_sub3(*X, *stap);
	su2_update_heat(X);

	//compute W13(2x2)
	su2_mat_mul_dag(W13, *X, A);

	//compute U' = stap = W13 * U2
	create_su3_sub3(W13_3x3, W13);
	su3_mat_mul(*stap, W13_3x3, U2);
	/*logging("\nnew link U'[%i][%i] = \n", n, dir);
	double *r;
	   r = (double*)(stap);
	   for(unsigned int i=0; i<18; i=i+2)
	   {
	       logging("(\%f,%f)", r[i], r[i+1]);
	       if(i==4 || i==10 || i==16)
	    	   logging("\n");
	       else
	    	   logging("\t");
	   }*/

	free(X);
}
#endif

