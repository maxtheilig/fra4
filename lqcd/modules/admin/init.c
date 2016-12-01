
/*******************************************************************************
*
* File init.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to initialise the fields and the error checks.
* 
* Externally accessible functions:
* 
* void init_program(int itype)
*      Initialises some smaller parts of the program, setting constants etc..
* 
* void neib_init(void)
*      Initialises the neib arrays.
* 
*******************************************************************************/

#define INIT_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"misc.h"
#include"headers.h"
#include"modules.h"

#define DEBUG_INIT 0



void init_program(int itype)
{
   int processes;
   
   error((runp.nc<=0)||(runp.dc<0)||(runp.tc<0)||
         ((((runp.nc-runp.tc)%runp.wc)!=0)&&(runp.wc>0)),
         "init_program [init.c]","Error in runparameters");
   
   PI=2.0*asin(1.0);
   
   init_twopi_start();
}



void neib_init(void)
{
   int ii,jj,kk,ll,ir,inei;
   int ix[DIM],ixn[DIM];
   int io[DIM];
   int ileng[DIM],iarea[DIM];
   
   checkpoint("neib_init -- in");
   
#if(DIM==2)
   io[0]=LENGS1;
   io[1]=1;
   ileng[0]=LENGT;
   ileng[1]=LENGS1;
   iarea[0]=TAREA;
   iarea[1]=LENGT;
#elif(DIM==3)
   io[0]=LENGS1*LENGS2;
   io[1]=LENGS2;
   io[2]=1;
   ileng[0]=LENGT;
   ileng[1]=LENGS1;
   ileng[2]=LENGS2;
   iarea[0]=TAREA;
   iarea[1]=LENGT*LENGS2;
   iarea[2]=LENGT*LENGS1;
#elif(DIM==4)
   io[0]=LENGS1*LENGS2*LENGS3;
   io[1]=LENGS2*LENGS3;
   io[2]=LENGS3;
   io[3]=1;
   ileng[0]=LENGT;
   ileng[1]=LENGS1;
   ileng[2]=LENGS2;
   ileng[3]=LENGS3;
   iarea[0]=TAREA;
   iarea[1]=LENGT*LENGS2*LENGS3;
   iarea[2]=LENGT*LENGS1*LENGS2;
   iarea[3]=LENGT*LENGS1*LENGS2;
#endif
   
   for(ii=0;ii<DIM;ii++)
   {
      lat.llength[ii]=ileng[ii];
      lat.barea[ii]=iarea[ii];
      lat.io[ii]=io[ii];
   }
   
   
   for(ii=0;ii<VOL;ii++)
   {
      for(jj=0;jj<DIM;jj++)
      {
	 neib[ii][jj]=-1;
	 neib[ii][jj+DIM]=-1;
      }
   }
   for(ii=0;ii<VOL;ii++)
   {
      for(jj=0,ir=ii;jj<DIM;jj++)
      {
         ix[jj]=(int)(ir/io[jj]);
         ir=ir%io[jj];
      }
      for(jj=0;jj<DIM;jj++)
      {
	 for(kk=0,ir=1;kk<DIM;kk++)
	 {
	    ixn[kk]=ix[kk];
	 }
	 
	 // Positive direction
	 ixn[jj]=ix[jj]+1;
	 if(ixn[jj]==ileng[jj])
	 {
	    ixn[jj]=0;
	    for(kk=0,neib[ii][jj]=0;kk<DIM;kk++)
	       neib[ii][jj]+=ixn[kk]*io[kk];
	 }
	 else
	 {
	    for(kk=0,neib[ii][jj]=0;kk<DIM;kk++)
	       neib[ii][jj]+=ixn[kk]*io[kk];
	 }
	 
	 // Negative direction
	 ixn[jj]=ix[jj]-1;
	 if(ixn[jj]<0)
	 {
	    ixn[jj]=ileng[jj]-1;
	    for(kk=0,neib[ii][jj+DIM]=0;kk<DIM;kk++)
	       neib[ii][jj+DIM]+=ixn[kk]*io[kk];
	 }
	 else
	 {
	    for(kk=0,neib[ii][jj+DIM]=0;kk<DIM;kk++)
	       neib[ii][jj+DIM]+=ixn[kk]*io[kk];
	 }
      }
   }
   
#if(DEBUG_INIT==1)
   logging("\nix");
   for(jj=0;jj<DIM;jj++)
      logging("\tn%d\tn%d",jj,jj+DIM);
   for(ii=0;ii<VOL;ii++)
   {
      logging("\n%d",ii);
      for(jj=0;jj<DIM;jj++)
	 logging("\t%d\t%d",neib[ii][jj],neib[ii][jj+DIM]);
   }
   logging("\n");
#endif
   checkpoint("neib_init -- out");
   return;
}
