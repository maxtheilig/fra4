
/*******************************************************************************
*
* File qcd.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program for
*
* Syntax: qcd -i [input-filename] -c [configuration file]
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include<stdio.h>
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
   update(1, 10, 0, 0);
   /*read_config("test_config");

   double *r;
   r = (double*)(pu[0][0]);
   for(unsigned int i=0; i<18; i=i+2)
   {
       logging("(\%f,%f)", r[i], r[i+1]);
       if(i==4 || i==10 || i==16)
    	   logging("\n");
       else
    	   logging("\t");
   }

   double x = plaquette();
   logging("plaquette = %f\n", x);
   double y = gauge_action();
   logging("S_G = %f\n", y);*/

   error_checks(1,0);
   print_info(seed,rconf,cnfg_file);
   
   finish_gauge();
   
   return 0;
}
