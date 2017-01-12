
/*******************************************************************************
*
* File spin_alg.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routine needed for spinor algebra.
* 
* Externally accessible functions:
* 
* void set_spinor_zero(sun_wferm *f)
*      Sets the spinor pointed to by f to 0.
* 
* void spin_copy(sun_wferm *f2,sun_wferm *f1)
*      Copies the spinor f1 into f2.
* 
* void spin_rmul_add(sun_wferm *f3,sun_wferm *f2,double a,sun_wferm *f1)
*      Computes f3=f2+c*f1 where f1,2,3 are the spinors pointed to by the
*      associated pointers and c is a real number.
* 
* double square_norm(sun_wferm *f)
*        Returns |f|^2 normalised to 1 for a unit spinor (i.e. the volume
*        factor 4*SUN*VOL has been divided out).
* 
* double spin_full_sum(sun_wferm *f)
*        Returns the sum of all components (including imaginary ones), again
*        normalised to 1 for a unit spinor.
* 
* complex scalar_prod(sun_wferm *f2,sun_wferm *f1)
*         Returns the scalar product (f2,f1) normalised to 1 when f2 and f1
*         are unit spinors.
* 
*******************************************************************************/

#define SPIN_ALG_C

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"ranlxd.h"
#include"headers.h"
#include"modules.h"

void set_spinor_ascending(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    int x = 1;
    for(;s<sf;s+=1)
    {
    	*s=x;
    	if(x==24)
    		x = 0;
    	++x;
    }
}

void set_spinor_one(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    int x = 0;
    for(;s<sf;s+=1, ++x)
    {
    	if(x%2==0)
    		*s=1.;
    	else
    		*s=0.;
    }
}

void set_spinor_zero(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    for(;s<sf;s+=1)
    {
       *s=0.;
    }
}



void spin_copy(sun_wferm *f2,sun_wferm *f1)
{
   double *s1,*s2,*sf;
   
   s1=(double*)(f1);
   s2=(double*)(f2);
   sf=s1+8*SUN*VOL;
   for(;s1<sf;s1+=1,s2+=1)
   {
      *s2=*s1;
   }
}



void spin_rmul_add(sun_wferm *f3,sun_wferm *f2,double a,sun_wferm *f1)
{
   double *s1,*s2,*s3,*sf;
   
   s1=(double*)(f1);
   s2=(double*)(f2);
   s3=(double*)(f3);
   sf=s1+8*SUN*VOL;
   for(;s1<sf;s1+=1,s2+=1,s3+=1)
   {
      *s3=*s2+a*(*s1);
   }
}



double square_norm(sun_wferm *f)
{
   double norm;
   complex *s,*sf;
   
   s=(complex*)(f);
   sf=s+4*SUN*VOL;
   for(norm=0.;s<sf;s+=1)
   {
      norm+=(*s).re*(*s).re+(*s).im*(*s).im;
   }
   
   return norm/(double)(4*SUN*VOL);
}



double spin_full_sum(sun_wferm *f)
{
   double sum;
   complex *s,*sf;
   
   s=(complex*)(f);
   sf=s+4*SUN*VOL;
   for(sum=0.;s<sf;s+=1)
   {
      sum+=(*s).re;
      sum+=(*s).im;
   }
   
   return sum/(double)(4*SUN*VOL);
}



complex scalar_prod(sun_wferm *f2,sun_wferm *f1)
{
   complex z,prod;
   complex *s1,*s2,*sf;
   
   s1=(complex*)(f1);
   s2=(complex*)(f2);
   sf=s1+4*SUN*VOL;
   prod.re=0.;
   prod.im=0.;
   for(;s1<sf;s1+=1,s2+=1)
   {
      compl_mult_sn(z,*s2,*s1);
      compl_selfadd(prod,z);
   }
   compl_realdiv_single(prod,(double)(4*SUN*VOL));
   
   return prod;
}
