
/*******************************************************************************
*
* File meas.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program for
*
* Syntax: meas -i [input-filename] -c [configuration file]
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
   init_gauge(0);
   //read_config("test_config");
   //specify run_inlmeas parameters
   inlmeas.ts = 4; inlmeas.tf = 7; inlmeas.dt = 4; inlmeas.rs = 4; inlmeas.rf = 7; inlmeas.dr = 4;
   inlmeas.ismt = 1; inlmeas.isms = 20; inlmeas.smpart = 0.5; inlmeas.smpars = 0.5;
   read_config("configs/24x24x24_SU2_b5.000000_id1_n1");
   //double w[1000];
   double *w;
   w = malloc(10*sizeof(double));
   meas_wilson(w);
   logging("W(4,4) = %f\n", w[0]);

   error_checks(1,0);
   print_info(seed,rconf,cnfg_file);

   finish_gauge();

   return 0;
}
