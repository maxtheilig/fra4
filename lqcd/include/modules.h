
/*******************************************************************************
*
* File modules.h
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Provides the prototypes for the global functions from all files.
*
*******************************************************************************/

#include<stddef.h>

#ifndef GAUGE_H
 #include"gauge.h"
#endif

#ifndef HEADERS_H
 #include"headers.h"
#endif

/* IO-archive functions */

#ifndef INP_IO_C
extern void read_inp(int*,char*,int*,char*,int,char *argv[]);
extern void setup_files(int,char*);
extern void print_info(int,int,char*);
#endif

#ifndef IO_UTILS_C
extern void logging(char*,...);
extern void error(int,char*,char*,...);
extern void checkpoint(char *format,...);
#endif

#ifndef CONFIG_IO_C
extern void write_config(char*);
extern void read_config(char*);
#endif

/* Initialisation */

#ifndef INIT_C
extern void init_program(int);
extern void init_gauge(int);
extern void finish_gauge(void);
extern void neib_init(void);
extern void allocate_gauge(sun_mat *u[VOL][DIM]);
extern void free_gauge(sun_mat *u[VOL][DIM]);
#endif

#ifndef RANDOM_SU3_C
extern void init_twopi_start(void);
extern int error_check_twopi_start(void);
extern void gauss(double*,int);
extern void random_su3_vector(su3vec*);
extern void random_su2(su2mat*);
extern void random_su3(su3mat*);
#endif

/* utils functions */

#ifndef UTILS_C
extern double get_time(void);
extern int custom_isnan(double);
extern int custom_isinf(double);
#endif

#ifndef ERROR_CHECKS_C
extern void set_warning(void);
extern void error_checks(int,int);
#endif

/* maths */

#ifndef SUN_VFUNC_C
extern void project_to_su3(su3mat *u);
extern void project_to_su2(su2mat *u);
extern void project_gfield_to_sun(sun_mat *u[VOL][DIM]);
#endif

/* observables */

#ifndef PLAQUETTE_C
extern double plaquette(void);
extern double gauge_action(void);
#endif
