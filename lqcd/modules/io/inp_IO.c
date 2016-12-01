
/*******************************************************************************
*
* File inp_IO.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines for the reading of the infile and the setup of the
* files, as well as some startup writeouts.
* 
* Externally accessile functions:
* 
* void read_inp(int *iseed,char *dir,int *rconf,
*               char *cfile,int argc,char *argv[])
*      Reads the input file and stores the associated parameters in the
*      struct 'run_para runp' (defined in headers.h). Some infos are
*      also pased back to the calling routine (typically the main) via
*      the arguments of the function.
* 
*      For parameters in run_para see the documentation.
*      iseed : seed for random number generator
*      dir: directory for output
*      rconf: Start from a given configuration?
*      cfile: Filename of that configuration
* 
* void setup_files(int id,char *dir)
*      Setup of the main files for output.
* 
*      id: run id
*      dir: output diectory
* 
* void print_info(int iseed,int rconf,char *cfile)
*      Printing the startup information into the logfile. for identification
*      of run parameters etc..
* 
*      iseed: seed of randum number generator
*      rconf,cfile: as above
* 
*******************************************************************************/

#define INP_IO_C

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"headers.h"
#include"modules.h"



void read_inp(int *iseed,char *dir,int *rconf,char *cfile,int argc,char *argv[])
{
   FILE *inf=NULL,*ftest=NULL;
   int i,ii,ifail;
   int id,nconf,nstep,ntherm,wcnfg;
   int ipl,ipc,irs,irf,ird;
   char itype[NAME_SIZE];
   double beta,eps;
   
   ftest=fopen("SETUP_ERRORS","r");
   if(ftest!=NULL)
   {
      fclose(ftest);
      remove("SETUP_ERRORS");
   }
   sprintf(LOG_FILE,"SETUP_ERRORS");
   
   for(i=1,ii=0;i<argc;i++)
      if (strcmp(argv[i],"-c")==0)
         ii=i+1;
   if(ii!=0)
   {
      strcpy(cfile,argv[ii]);
      inf=fopen(cfile,"r");
      error(inf==NULL,"read_inp [inp_IO.c]",
            "Incorrect config file!\n"
            "Syntax: IND-G_pg -i <input file> [-c <input config>]");
      fclose(inf);
      inf=NULL;
      *rconf=1;
   }
   else
      *rconf=0;
   
   for(i=1,ii=0;i<argc;i++)
      if (strcmp(argv[i],"-i")==0)
         ii=i+1;
   error(ii==0,"read_inp [inp_IO.c]",
         "Syntax: IND-G_pg -i <input file> [-c <input config>]");
   inf=fopen(argv[ii],"r");
   error(inf==NULL,"read_inp [inp_IO.c]",
         "Unable to open input file!");
   
   ifail=fscanf(inf,"id      %d\n",&id);
   ifail=fscanf(inf,"odir    %s\n",dir);
   ifail=fscanf(inf,"beta    %lf\n",&beta);
   ifail=fscanf(inf,"epsilon %lf\n",&eps);
   ifail=fscanf(inf,"nconf   %d\n",&nconf);
   ifail=fscanf(inf,"ntherm  %d\n",&ntherm);
   ifail=fscanf(inf,"nstep   %d\n",&nstep);
   ifail=fscanf(inf,"wconfig %d\n",&wcnfg);
   ifail=fscanf(inf,"seed    %d\n",iseed);
   
   fclose(inf);
   
   runp.id=id;
   runp.beta=beta;
   runp.eps=eps;
   runp.nc=nconf;
   runp.dc=nstep;
   runp.tc=ntherm;
   runp.wc=wcnfg;
}



void setup_files(int id,char *dir)
{
   FILE *fset=NULL;
   char base[NAME_SIZE],check[NAME_SIZE];
   
   if(DIM==2)
      sprintf(base,"%dx%d_SU%d_b%.6f_id%d",LENGT,LENGS1,
	      SUN,runp.beta,id);
   else if(DIM==3)
      sprintf(base,"%dx%dx%d_SU%d_b%.6f_id%d",LENGT,LENGS1,LENGS2,
	      SUN,runp.beta,id);
   else
      sprintf(base,"%dx%dx%dx%d_SU%d_b%.6f_id%d",LENGT,LENGS1,LENGS2,LENGS3,
	      SUN,runp.beta,id);
   
   sprintf(check,"%s/%s.log",dir,base);
   fset=fopen(check,"r");
   error(fset!=NULL,"setup_files [inp_IO.c]",
         "Try to overwrite existing .log file!");
   sprintf(check,"%s/%s.out",dir,base);
   fset=fopen(check,"r");
   error(fset!=NULL,"setup_files [inp_IO.c]",
         "Try to overwrite existing .out file!");
   
   fset=fopen(LOG_FILE,"r");
   if(fset!=NULL)
   {
      printf("SETUP_ERRORS file was not empty!!!\n");
      printf("Aborting!\n");
      exit(-1);
   }
   
   sprintf(LOG_FILE,"%s/%s.log",dir,base);
   sprintf(OUT_FILE,"%s/%s.out",dir,base);
   sprintf(CNFG_FILE,"%s/%s",dir,base);
   
   fset=fopen(LOG_FILE,"w");
   error(fset==NULL,"setup_files [inp_IO.c]",
         "Unable to create .log file!");
   fclose(fset);
   fset=fopen(OUT_FILE,"w");
   error(fset==NULL,"setup_files [inp_IO.c]",
         "Unable to create .out file!");
   fclose(fset);
}



void print_info(int iseed,int rconf,char *cfile)
{
   FILE *flog=NULL;
   int ii;
   
   flog=fopen(LOG_FILE,"ab");
   error(flog==NULL,"print_info_heat [inp_IO.c]",
         "Unable to open .log file!");
   
   fprintf(flog,"*******************************\n");
   fprintf(flog,"LATTICE QCD CODE\n");
   fprintf(flog,"*******************************\n");
   fprintf(flog,"Global parameters:\n");
   fprintf(flog,"*******************************\n");
#if(DIM==2)
   fprintf(flog,"Using a two-dimensional lattice:\n");
   fprintf(flog,"latt  :\t%dx%d\n",LENGT,LENGS1);
#elif(DIM==3)
   fprintf(flog,"Using a three-dimensional lattice:\n");
   fprintf(flog,"latt  :\t%dx%dx%d\n",LENGT,LENGS1,LENGS2);
#else
   fprintf(flog,"Using a four -dimensional lattice:\n");
   fprintf(flog,"latt  :\t%dx%dx%dx%d\n",LENGT,LENGS1,LENGS2,LENGS3);
#endif
   fprintf(flog,"beta :\t%.6f\n",runp.beta);
   fprintf(flog,"eps   :\t%.6f\n",runp.eps);
   fprintf(flog,"nconf :\t%d\n",runp.nc);
   fprintf(flog,"dsw   :\t%d\n",runp.dc);
   fprintf(flog,"therm :\t%d\n",runp.tc);
   fprintf(flog,"nsave :\t%d\n",runp.wc);
   fprintf(flog,"seed  :\t%d\n",iseed);
   fprintf(flog,"*******************************\n");
   fprintf(flog,"Inline measurements:\n");
   fprintf(flog,"None at the moment!\n");
   fprintf(flog,"*******************************\n");
   fprintf(flog,"Start from ??? configuration:\n");
   fprintf(flog,"*******************************\n");
   
   fclose(flog);
}
