
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
#include <time.h>


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
   double x;
   int stype = 0;
   int utype = 0; //0=metro, 1=heat
   FILE *write;
   write = fopen("plaq_metro_su2.txt", "w");
   for(unsigned int i=0; i<runp.tc; ++i)
   {
	   update(1, 1, utype, stype);
	   x = plaquette();
	   logging("plaq = %f\n", x);
	   fprintf(write, "%f\n", x);
   }
   //clock_t prgstart, prgende; //Definierung der Variablen
   //prgstart=clock(); //CPU-Zeit zu Beginn des Programmes
   for(unsigned int i=0; i<runp.nc; ++i)
   {
	   /*for(unsigned int k=0; k<runp.wc; ++k)
		   update(1, 1, utype, stype);*/
	   update(1, 1, utype, stype);
	   x = plaquette();
	   logging("plaq = %f\n", x);
	   fprintf(write, "%f\n", x);
   }
   fclose(write);
   //prgende=clock();//CPU-Zeit am Ende des Programmes
   //printf("Laufzeit %.8f Sekunden\n",(float)(prgende-prgstart) / CLOCKS_PER_SEC);
   //read_config("test_config");

   /*double *r;
   r = (double*)(pu[0][0]);
   for(unsigned int i=0; i<18; i=i+2)
   {
       logging("(\%f,%f)", r[i], r[i+1]);
       if(i==4 || i==10 || i==16)
    	   logging("\n");
       else
    	   logging("\t");
   }*/

   error_checks(1,0);
   print_info(seed,rconf,cnfg_file);
   
   finish_gauge();
   
   return 0;
}
