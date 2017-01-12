
/*******************************************************************************
*
* File error_checks.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes some error checks while running the program.
* 
* Externally accessible functions:
* 
* void set_warning(void)
*      Increases the internal warning counter warn_log.
* 
* void error_checks(int safe)
*      Provides internal consistency checks for the state of the program.
*      If safe!=0 the status of the warning counter is checked and the
*      program aborted in case of any present warnings.
* 
*******************************************************************************/

#define ERROR_CHECKS_C

#include<stdlib.h>
#include<stdio.h>
#include<float.h>
#include<math.h>
#include"headers.h"
#include"modules.h"

static int warn_log=0;



void set_warning(void)
{
   warn_log++;
}



void error_checks(int safe)
{
   int iam_rank;
   
   error(fabs(PI-2.0*asin(1.0))>1.E-15,"error_checks",
	 "PI not initialised!");
   error(error_check_twopi_start(),"error_checks",
	 "twopi not initialised in random_su3!");
   
   if(safe==1)
   {
      error(warn_log!=0,"error_checks",
	    "There were warnings!\nProgram safely aborted!");
   }
}
