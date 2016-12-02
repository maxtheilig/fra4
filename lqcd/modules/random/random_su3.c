
/*******************************************************************************
*
* File random_su3.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to generate a random SU(3) matrix and bosons.
* 
* The routines included are similar to those used by Martin Luescher in
* the DD-HMC code.
* 
* Externally accessible functions:
* 
* void init_twopi_start(void)
*      Initialises the local twopi variable.
* 
* int error_check_twopi_start(void)
*     Checks whether twopi has been initialised.
* 
* void gauss(double *rd,int n)
*      Returns n Gaussian distributed (e^{-r^2}) random numbers in the
*      array rd.
* 
* void random_su3_vector(su3vec *v)
*      Returns a random su3vec v with non-zero norm, including Gaussian
*      distributed random numbers.
*
* void random_su2(su2mat *u)
*      Generates a random su2mat using a random su3vec.
*
* void random_su3(su3mat *u)
*      Generates a random su3mat using three random su3vec.
*
*******************************************************************************/

#define RANDOM_SU3_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include"headers.h"
#include"modules.h"
#include"ranlxd.h"

static int init_twopi=0;
static double twopi;



void init_twopi_start(void)
{
   twopi=2.*PI;
   init_twopi=1;
}



int error_check_twopi_start(void)
{
   return (init_twopi==0);
}



void gauss(double *rd,int n)
{
   int k;
   double ud[2];
   double x1,x2,rho,y1,y2;

   for(k=0;k<n;)
   {
      ranlxd(ud,2);
      x1=ud[0];
      x2=ud[1];
      rho=sqrt(-log(1.0-x1));
      x2*=twopi;
      y1=rho*sin(x2);
      y2=rho*cos(x2);

      rd[k++]=y1;
      if(k<n)
         rd[k++]=y2;
   }
}



void random_su3_vector(su3vec *v)
{
   int i;
   double *r,norm,fact;

   r=(double*)(v);
   norm=0.0;

   while (norm<=DBL_EPSILON)
   {
      gauss(r,6);
      norm=0.0;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=sqrt(norm);
   }

   fact=1.0/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;
}

void random_su2(su2mat *u)
{
	logging("random_su2 not implemented yet");
/*	su3vec *v;
	random_su3_vector(v);
	u->c0 = sqrt(1 - v->c1.re*v->c1.re - v->c2.re*v->c2.re - v->c3.re*v->c3.re);
	u->c1 = v->c1.re;
	u->c2 = v->c2.re;
	u->c3 = v->c3.re;*/
}

void random_su3(su3mat *u)
{
	su3vec *v1, *v2, *v3;
	double *r, norm, fact;

	v1 = (su3vec*)(u);
	v2 = v1 + 1;
	v3 = v1 + 2;

	r=(double*)(v3);

	norm=0.0;
	while (norm<=DBL_EPSILON)
	{
		random_su3_vector(v1);
		random_su3_vector(v2);
		su3vec_cross_prod(*v3, *v1, *v2);
		norm=0.0;
		for (unsigned int i=0; i<6; i++)
			norm+=r[i]*r[i];
		norm=sqrt(norm);
	}

    fact=1.0/norm;
    for (unsigned int i=0; i<6; i++)
       r[i]*=fact;

    su3vec_cross_prod(*v2, *v3, *v1);
}
