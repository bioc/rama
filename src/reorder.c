#include "util.h" 

void reorder(double *data_vec, int *n1, int *n2, double *all_data, int *m1, int *m2)
{
  int i,j,x,y;
  double **data;
  
  data=dmatrix(*n1, *n2);  
  vec_mat(data_vec,n1,n2,data);
  
  for(i=0;i<*n1;i++)
    for(j=0;j<*n2;j++)
      {
	/** Coordinates in the right order **/
	x=data[i][0];
	y=data[i][1];
	
	/** Recopy the matrix at the right place **/
	all_data[(x**m2+y)**n2+j]=data[i][j];
      }
  
  free_dmatrix(data, *n1);
  
}
 
