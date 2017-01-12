
/*******************************************************************************
*
* File update.c
*
* Copyright (C) 2016 Bastian Brandt, Max Theilig
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Externally accessible functions:
*
*******************************************************************************/

#define UPDATE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include"headers.h"
#include"modules.h"
#include "ranlxd.h"

void update(int iup, int nup, int utype, int stype)
{
	double acc = 0.0;
	double r[2];
	int x,y;
	sun_mat *stap;
	stap = malloc(sizeof(sun_mat));
	if(utype == 0){
		if(stype == 0){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int n=0; n<VOL; ++n)
				{
					for(unsigned int dir=0; dir<DIM; ++dir)
					{
							acc += local_metr(n, dir, 1);
					}
				}
			}
		}
		if(stype == 1){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int i=0; i<(VOL*DIM); ++i)
				{
					ranlxd(r, 2);
					x = (int) (r[0] * VOL);
					if(x == VOL) { x = x - 1; }
					y = (int) (r[1] * DIM);
					acc += local_metr(x, y ,1);
				}
			}
		}
		acc /= VOL*DIM*nup;
		//logging("Acceptance Propability = %f\n", acc);
	}

	if(utype == 1){
#if(SUN==2)
		if(stype == 0){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int n=0; n<VOL; ++n)
				{
					for(unsigned int dir=0; dir<DIM; ++dir)
					{
						staples(n, dir, stap);
						su2_update_heat(stap);
						*pu[n][dir] = *stap;
					}
				}
			}
		}
		if(stype == 1){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int i=0; i<(VOL*DIM); ++i)
				{
					ranlxd(r, 2);
					x = (int) (r[0] * VOL);
					if(x == VOL) { x = x - 1; }
					y = (int) (r[1] * DIM);
					staples(x, y, stap);
					su2_update_heat(stap);
					*pu[x][y] = *stap;
				}
			}
		}
		project_gfield_to_sun(pu);
		for(unsigned int j=0; j<runp.no; ++j)
		{
			for(unsigned int n=0; n<VOL; ++n)
			{
				for(unsigned int dir=0; dir<DIM; ++dir)
				{
					staples(n, dir, stap);
					su2_relax(n, dir, stap);
					*pu[n][dir] = *stap;
				}
			}
		}
#elif(SUN==3)
		if(stype == 0){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int n=0; n<VOL; ++n)
				{
					for(unsigned int dir=0; dir<DIM; ++dir)
					{
						staples(n, dir, stap);
						su3_update_heat(n, dir, stap);
						//project_to_su3(stap);
						*pu[n][dir] = *stap;
					}
				}
			}
		}
		if(stype == 1){
			for(unsigned int k=0; k<nup; ++k)
			{
				for(unsigned int i=0; i<(VOL*DIM); ++i)
				{
					ranlxd(r, 2);
					x = (int) (r[0] * VOL);
					if(x == VOL) { x = x - 1; }
					y = (int) (r[1] * DIM);
					staples(x, y, stap);
					su3_update_heat(x, y, stap);
					//project_to_su3(stap);
					*pu[x][y] = *stap;
				}
			}
		}
#endif
	}
	project_gfield_to_sun(pu);
	free(stap);
}
