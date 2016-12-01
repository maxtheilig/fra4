
/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes some additional routines.
* 
* Externally accessible functions:
* 
* double get_time(void)
*        Returns the current execution time after the first call of get_time.
*        I.e. the first call return 0!
* 
* int custom_isnan(double var)
*     Returns true when var is not-a-number.
* 
* int custom_isinf(double x)
*     Returns true if x is infinite.
* 
*******************************************************************************/

#define UTILS_C

#include<stdlib.h>
#include<stdio.h>
#include"headers.h"
#include"modules.h"
#include<time.h>

static int time_flag=0;
static time_t ref_time;



double get_time(void)
{
   if(time_flag==0)
   {
      time(&ref_time);
      time_flag=1;
      return 0.;
   }
   else
   {
      time_t wt;
      time(&wt);
      return difftime(wt,ref_time);
   }
}



int custom_isnan(double var)
{
    volatile double d = var;
    return d != d;
}



int custom_isinf(double x)
{
   volatile double temp = x;
   if ((temp == x) && ((temp - x) != 0.0))
      return 1;
   else return 0;
}
