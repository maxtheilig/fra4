
/*******************************************************************************
*
* File headers.h
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* For description of the parameters see the file README.main in ../main
*
*******************************************************************************/

#define HEADERS_H

#define DIM 4
#define SUN 3

#define LENGT 8
#define LENGS1 8
#define LENGS2 8
#define LENGS3 8

#define NAME_SIZE 128
#define DEBUG 0

/* SIM_TYPE==0 : Metropolis */
#define SIM_TYPE 0

/******************************************************************************
* Checks
*******************************************************************************/

#if((DIM<2)||(DIM>4)||(SUN<2)||(SUN>3))
#error : DIM or SUN not suitable
#endif

#if((LENGT<2)||(LENGS1<4)||(LENGS2<4)||(LENGS3<4))
#error : The lattice is too small
#endif

#if(((LENGT%2)!=0)||((LENGS1%2)!=0)||((LENGS2%2)!=0)||((LENGS3%2)!=0))
#error : The lattice has to consist of an even number of points
#endif

#if((DEBUG!=0)&&(DEBUG!=1))
#error : DEBUG must be set to 0 or 1
#endif

#if(NAME_SIZE<128)
#error : NAME_SIZE should be bigger or equal to 128
#endif

/******************************************************************************
* Derived and additional definitions
*******************************************************************************/

#if(DIM==2)
#define VOL (LENGT*LENGS1)
#define TAREA LENGS1
#define NPLAQ 1
#elif(DIM==3)
#define VOL (LENGT*LENGS1*LENGS2)
#define TAREA (LENGS1*LENGS2)
#define NPLAQ 3
#elif(DIM==4)
#define VOL (LENGT*LENGS1*LENGS2*LENGS3)
#define TAREA (LENGS1*LENGS2*LENGS3)
#define NPLAQ 6
#endif

/******************************************************************************
* Definition of gauge functions:
*******************************************************************************/

/* Gauge macros */
#ifndef GAUGE_H
 #include"gauge.h"
#endif

#if(SUN==2)
 #define SUNVOL 4
 #define ALGVOL 3
 #define sun_mat su2mat
 #define sun_alg su2alg
 #define sun_dag su2_dag
 #define sun_star su2_star
 #define sun_star_2g su2_star_2g
 #define sun_transp su2_transp
 #define sun_unit su2_unit
 #define sun_zero su2_zero
 #define sun_mul su2_mat_mul
 #define sun_mul_dag su2_mat_mul_dag
 #define sun_add su2_add
 #define sun_self_add su2_self_add
 #define sun_sub su2_sub
 #define sun_self_sub su2_self_sub
 #define sun_trace su2_trace
 #define sun_trace_cl su2_trace_cl
 #define mksun_nxn mksu2_2x2
 #define mksun_nxn_dag mksu2_2x2_dag
 #define sun_dble_mul su2_dble_mul
 #define sun_dble_div su2_dble_div
 #define mk_sun_dble_array mk_su2_dble_array
 #define mk_dble_array_sun mk_dble_array_su2
 #define mk_dble_array_sun_alg mk_dble_array_su2_alg
 #define mk_sun_alg_dble_array mk_su2_alg_dble_array
 #define set_alg_zero set_alg_zero_su2
 #define PROJ_FREQ 10
 #define random_sun random_su2
#elif(SUN==3)
 #define SUNVOL 18
 #define ALGVOL 8
 #define sun_mat su3mat
 #define sun_alg su3alg
 #define sun_dag su3_dag
 #define sun_star su3_star
 #define sun_star_2g su3_star_2g
 #define sun_transp su3_transp
 #define sun_unit su3_unit
 #define sun_zero su3_zero
 #define sun_mul su3_mat_mul
 #define sun_mul_dag su3_mat_mul_dag
 #define sun_add su3_add
 #define sun_self_add su3_self_add
 #define sun_sub su3_sub
 #define sun_self_sub su3_self_sub
 #define sun_trace su3_trace_re
 #define sun_trace_cl su3_trace
 #define mksun_nxn mksu3_3x3
 #define mksun_nxn_dag mksu3_3x3_dag
 #define sun_dble_mul su3_dble_mul
 #define sun_dble_div su3_dble_div
 #define mk_sun_dble_array mk_su3_dble_array
 #define mk_dble_array_sun mk_dble_array_su3
 #define mk_dble_array_sun_alg mk_dble_array_su3_alg
 #define mk_sun_alg_dble_array mk_su3_alg_dble_array
 #define set_alg_zero set_alg_zero_su3
 #define random_sun random_su3
#endif

/******************************************************************************
* Global variables and structures:
*******************************************************************************/

typedef struct
{
   int id,nc,tc,dc,niter,nw,wc;
   double m,beta,tau,eps;
} run_para;

typedef struct
{
   int imeas;
} run_inlmeas;

typedef struct
{
   int llength[DIM];
   int barea[DIM];
   int io[DIM];
} lat_parms;

#ifdef MAIN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN run_para runp;
EXTERN run_inlmeas inlmeas;
EXTERN lat_parms lat;

EXTERN sun_mat *pu[VOL][DIM];

EXTERN int neib[VOL][2*DIM];

EXTERN char LOG_FILE[NAME_SIZE];
EXTERN char OUT_FILE[NAME_SIZE];
EXTERN char CNFG_FILE[NAME_SIZE];

EXTERN double PI;

#undef EXTERN
