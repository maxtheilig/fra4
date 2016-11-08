
/*******************************************************************************
*
* File autocorr.c
*
* Copyright (C) 2016 Sebastian Schmalzbauer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Statistical analysis of data from a file
*
*  compile like:
*  cc -O2 -DSSE autocorr.c -lm -o autocorr
*
*  run:
*  ./autocorr <data.txt>
*
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double average(double *array, int start, int end)
{
   int i;
   double average = 0.0;
   for(i = start; i < end; i++)
      average += array[i];
   return average/((double)(end - start));
}


double covariance_k(double *array, int start, int end, int k)
{
   int i;
   double av0 = average(array, start, end-k);
   double avk = average(array, start + k, end);
   double c = 0.0;
   for(i = 0; i < end - k; i++)
      c += array[i]*array[i + k];
   return c/((double)(end - start - k)) - av0*avk;
}


double tau_int(double *array, int start, int end)
{
   int k;
   double av = average(array, start, end);
   double var = covariance_k(array, start, end, 0);
   double t, c;
         
   if(var == 0.0)
      return 0.5;
   else
   {
      t = var/2.0;
      for(k=1; k < end; k++)
      {
         c = covariance_k(array, start, end, k);
         if(c<0.0)
            break;
         t += c;
      }      
      return t/var;
   }
}


double correlated_error(double *array, int start, int end)
{
   return sqrt(covariance_k(array, start, end, 0)/(double)(end - start - 1.0));
}


double uncorrelated_error(double *array, int N)
{
   double t = tau_int(array, 0, N);
   int rel = (int)ceil(20.0*t);
      
   if(rel>N/5 || rel == (int)ceil(NAN) || rel == (int)ceil(INFINITY))
   {
      printf("Data not equilibrated or need more measurements!\n");
      return INFINITY;
   }
   else
   {
      double uncorr_error = sqrt(2.0*t)*correlated_error(array, 0, N);
      printf("%f  +-  %f    with tau_int = %f\n", average(array, 0, N), uncorr_error, 2.0*t);
      return uncorr_error;
   }
}


int check_file(char *filename)
{
   double buff;
   int lines = 0;
      
   FILE *check;
   if (check = fopen(filename, "r"))
   {
      while (fscanf(check, "%lf\n", &buff) != EOF)
	 lines ++;
      fclose(check);
      return lines;
   }
   else
   {
      printf("Error: could not open %s!\n", filename);
      return -1;
   }
}


void read_data(char *filename, double *data)
{
   int i = 0;
   FILE *read = fopen(filename, "r");
   while (fscanf(read, "%lf\n", &data[i]) != EOF)
      i++;
   fclose(read);
}


int main(int argc,char *argv[])
{
   int i;
   if(argc != 2)
   {
      printf("Specify parameters: ./%s <data.txt>\n", argv[0]);
      exit(1);
   }
   
   int N = check_file(argv[1]);
   printf("file %s has %i lines\n", argv[1], N);
   
   double *data = malloc(N*sizeof(double));
   read_data(argv[1], data);
   
   uncorrelated_error(data, N); 
   exit(0);
}
