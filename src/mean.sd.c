#include "util.h"

void mean_sd(double **data, int *n1, int *n2, double *mean, double *sd);


/****************************************************************************************************/
/*                                                                   link_R_mean_sd                                                                               */
/*  Purpose:   Interface between mean_sd (C) and R                                                                                             */
/* Argument description:                                                                                                                                       */
/* Input:                                                                                                                                                                  */
/* data_vec: The data matrix in a vector form                                                                                                        */
/* n1: the number of rows of the matrix                                                                                                                */
/* n2: The number of columns of the matrix                                                                                                          */
/* Output:                                                                                                                                                                */
/* mean: The vector of means for all the genes                                                                                                       */
/* sd: The vector of sds for all the genes                                                                                                                 */
/*****************************************************************************************************/ 

void link_R_mean_sd(double *data_vec, int *n1, int *n2, double *mean, double *sd)
{
  
  double **data;
  
  /** Allocate the memory **/
  data=dmatrix(*n1,*n2);

  /** Reshape the data into matrices **/
  vec_mat(data_vec,n1,n2,data);

  /** Compute the means and the sds for every genes **/
  mean_sd(data, n1, n2, mean, sd);

  /** Free the tmp matrix **/
  free_dmatrix(data, *n1);
  
}

/* This function computes the means and standard deviations of all the genes (see description above) */ 
void mean_sd(double **data, int *n1, int *n2, double *mean, double *sd)
{
  int i;
  int no_use;

  for(i=0;i<*n1;i++)
    {
      mean[i]=mean_vec(data[i],n2);
      sd[i]=stdd(data[i], n2, &no_use); 
    }
}
