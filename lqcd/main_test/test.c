
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
   
   error_checks(1,0);
   print_info(seed,rconf,cnfg_file);
   
   
   
   return 0;
}
