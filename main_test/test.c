
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
   
   runp.m = 1.0;

   init_gauge(0);
   double *p; //ascending gaugefield
   for(unsigned int n=0; n<VOL; ++n){
	   for(unsigned int dir=0; dir<DIM; ++dir){
		  p = (double*)(pu[n][dir]);
		  for(unsigned int i=0; i<18; ++i)
		  {
			  p[i] = i+1;
		  }
	   }
   }

   sun_wferm *s[VOL], *r[VOL];
   alloc_wferm(s);
   alloc_wferm(r);
   set_spinor_one(*s);
   apply_dirac_wil(*r, *s);
   double x = square_norm(*r);
   //double x = spin_full_sum(*r);
   x *= 4*SUN;
   logging("x = %f\n", x);


   free_wferm(s);
   free_wferm(r);

   error_checks(1,0);
   print_info(seed,rconf,cnfg_file);
   
   return 0;
}
