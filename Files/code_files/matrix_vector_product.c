#include<stdio.h>

double matrix_vector_product(int n,int N_j,double ap[n], double aw[n], double as[n], double an[n], double ae[n], double d[n], double Ad[n])
{

	
		for(int i=0;i<n;i++)
		{
			//top row of A matrix 111111111111111111111111111111
			if(i==0)
			{
				Ad[i] = ap[i]*d[i] + an[i]*d[i+1] + ae[i]*d[i+N_j]; 
			} 
			
			//till start of west coefficient 22222222222222222222222222222222
			else if(i>0 && i<N_j)
			{
				Ad[i] = as[i]*d[i-1] + ap[i]*d[i] + an[i]*d[i+1] + ae[i]*d[i+N_j];
			}
			
			//from starting of west coefficient till end of east coefficient 3333333333333333333333333333
			else if(i>=N_j && i<n-N_j)
			{
				Ad[i] = aw[i]*d[i-N_j] + as[i]*d[i-1] + ap[i]*d[i] + an[i]*d[i+1] + ae[i]*d[i+N_j];
			}
			
			///from the end of east coefficient to the till end of north coefficient 444444444444444444444444444444444
			else if(i>=n-N_j && i<n-1)
			{
				Ad[i] = aw[i]*d[i-N_j] + as[i]*d[i-1] + ap[i]*d[i] + as[i]*d[i+1];
			}
			
			//last row of A matrix
			else if(i==n-1)
			{
				Ad[i] = aw[i]*d[i-N_j] + as[i]*d[i-1] + ap[i]*d[i];
			}
			
		}
		for(int i=0;i<n;i++)
		{
			//printf("Ad[%d]=%lf\n",i,Ad[i]);
		}
		
		
}
