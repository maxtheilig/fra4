
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
/* Observables */

#ifndef PLAQUETTE_C
extern double plaquette(void);
extern double gauge_action(void);
#endif

   /* New part */
#ifndef SMEARING_C
extern void smearing_APE_spatial(int,double,sun_mat *u[VOL][DIM]);
extern void smearing_APE_temporal(int,double,sun_mat *u[VOL][DIM]);
#endif

#ifndef WILSON_C
extern void meas_wilson(double*);
#endif
   /* */

/* Updates */

#ifndef EXP_FCT_C
extern void expx(double,sun_alg*,sun_mat*);
#endif

#ifndef METRO_C
extern void staples(int,int,sun_mat*);
extern double local_metr(int,int,int);
#endif

#ifndef UPDATE_C
extern void update(int,int,int,int);
#endif

#ifndef HEAT_C
extern void su2_update_heat(su2mat*);
extern void su3_update_heat(int,int,su3mat*);
extern void su2_relax(int,int,su2mat*);
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
extern void neib_init(void);
extern void init_gauge(int);
extern void finish_gauge(void);
extern void allocate_gauge(sun_mat *u[VOL][DIM]);
extern void free_gauge(sun_mat *u[VOL][DIM]);
extern void copy_gauge_field(sun_mat *u1[VOL][DIM],sun_mat *u2[VOL][DIM]);
extern void alloc_wferm(sun_wferm**);
extern void free_wferm(sun_wferm**);
#endif

#ifndef RANDOM_SU3_C
extern void init_twopi_start(void);
extern int error_check_twopi_start(void);
extern void gauss(double*,int);
extern void random_su3_vector(su3vec*);
extern void random_su3(su3mat*);
extern void random_su2(su2mat*);
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

/* dirac */
#ifndef DIRAC_WIL_C
extern void apply_dirac_wil(sun_wferm*, sun_wferm*);
#endif

#ifndef SPIN_ALG_C
extern void set_spinor_ascending(sun_wferm*);
extern void set_spinor_one(sun_wferm*);
extern void set_spinor_zero(sun_wferm*);
extern void spin_copy(sun_wferm* ,sun_wferm*);
extern void spin_rmul_add(sun_wferm* ,sun_wferm* ,double, sun_wferm*);
extern double square_norm(sun_wferm*);
extern double spin_full_sum(sun_wferm*);
extern complex scalar_prod(sun_wferm*, sun_wferm*);
#endif
