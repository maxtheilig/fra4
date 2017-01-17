
/*******************************************************************************
*
* File test.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program for
*
* Syntax: test -i <input-filename>
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include<stdio.h>
#include<stdlib.h>
#include"ranlxd.h"
#include"headers.h"
#include"modules.h"
#include"math.h"



void set_spinor_cold(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    for(;s<sf;s+=2)
    {
       *s=1.;
       *(s+1)=0.;
    }
}



void set_spinor_asc(sun_wferm *f)
{
    int ii;
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    for(;s<sf;)
       for(ii=0;ii<8*SUN;ii++,s++)
          *s=(double)(ii+1);
}



void set_unit_gauge(void)
{
   int n,dir;
   for(n=0;n<VOL;n++)
      for(dir=0;dir<DIM;dir++)
      {
	 sun_unit(*pu[n][dir]);
      }
}



void set_asc_gauge(void)
{
   int n,dir,ii;
   double *du;
   
   for(n=0;n<VOL;n++)
      for(dir=0;dir<DIM;dir++)
      {
	 du=(double*)(pu[n][dir]);
	 for(ii=0;ii<SUNVOL;ii++,du++)
	    *du=(double)(ii+1);
      }
}



int main(int argc,char *argv[])
{
   int seed,rconf;
   char out_dir[NAME_SIZE];
   char cnfg_file[NAME_SIZE];
   
   read_inp(&seed,out_dir,&rconf,cnfg_file,argc,argv);
   rlxd_init(1,seed);   
   init_program(1);
   
   setup_files(runp.id,out_dir);
   neib_init();
   init_gauge(1);
   
   if(rconf)
      read_config(cnfg_file);
   
   double norm,sum;
   sun_wferm *f1,*f2,*f3,*f4;
   alloc_wferm(&f1);
   alloc_wferm(&f2);
   alloc_wferm(&f3);
   alloc_wferm(&f4);
   
   /* Don't forget to comment/uncomment the antiperiodic boundary conditions in modules/dirac/dirac_wil.c  */
   
   runp.m=0.;
   logging("\nm=0:\n\n");
   
   set_spinor_cold(f1);
   set_unit_gauge();   
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_1, phi_1:   Psi^dag(D^dag D Psi)=%.0f\t\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_asc(f1);
   set_unit_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_1, phi_2:   Psi^dag(D^dag D Psi)=%.0f\t\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_cold(f1);
   set_asc_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_2, phi_1:   Psi^dag(D^dag D Psi)=%.0f\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_asc(f1);
   set_asc_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_2, phi_2:   Psi^dag(D^dag D Psi)=%.0f\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));


   runp.m=1.;
   logging("\nm=1:\n\n");

   set_spinor_cold(f1);
   set_unit_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_1, phi_1:   Psi^dag(D^dag D Psi)=%.0f\t\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_asc(f1);
   set_unit_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_1, phi_2:   Psi^dag(D^dag D Psi)=%.0f\t\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_cold(f1);
   set_asc_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_2, phi_1:   Psi^dag(D^dag D Psi)=%.0f\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));

   set_spinor_asc(f1);
   set_asc_gauge();
   apply_dirac_wil(f2,f1);
   apply_dirac_wil_conj(f3,f2);
   logging("U_2, phi_2:   Psi^dag(D^dag D Psi)=%.0f\t(D Psi)^dag (D Psi)=%.0f\n",
		   4*SUN*scalar_prod(f1,f3).re, 4*SUN*square_norm(f2));


   set_spinor_zero(f2);
   point_source(f1, 0, 0);

   logging("\n\n\nPoint Source Checks:\n");
   logging("square norm of point source: %2.2f\n", VOL*4*SUN*square_norm(f1));
   logging("value at n=0, d=0: %2.2f\n\n\n", f1[0].d1.c1.re );

   set_unit_gauge();
   runp.m = 2.;
   double eps = pow(10,-11), Xcond;
   int nmax = 1000;

   cg(f2, apply_dirac_wil, apply_dirac_wil_conj, f1, eps, nmax);

   apply_dirac_wil(f3, f2);

   sunwferm_sub(*f4, *f1, *f3);

   logging("square norm of error: %2.8f\n", square_norm(f4));
   if(square_norm(f4) < eps)
	   logging("error smaller than eps!\n");
   else
	   logging("error bigger than eps!\n");

   Xcond = scalar_prod(f1, f2).re;

   logging("value of X condensate: %2.5f\n", Xcond);

   free_wferm(&f1);
   free_wferm(&f2);
   free_wferm(&f3);
   free_wferm(&f4);
   
   return 0;
}
