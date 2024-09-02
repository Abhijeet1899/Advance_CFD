#include<stdio.h>

double residual(int n,int N_j,double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double x[n], double b[n])
{
	double h=1.0/n;
	
		for(int i=0;i<n;i++)
		{
			//top row of A matrix 111111111111111111111111111111
			if(i==0)
			{
				r[i] = b[i] - ap[i]*x[i] - an[i]*x[i+1] - ae[i]*x[i+N_j]; 
			} 
			
			//till start of west coefficient 22222222222222222222222222222222
			else if(i>0 && i<N_j)
			{
				r[i] = b[i] - as[i]*x[i-1] - ap[i]*x[i] - an[i]*x[i+1] - ae[i]*x[i+N_j];
			}
			
			//from starting of west coefficient till end of east coefficient 3333333333333333333333333333
			else if(i>=N_j && i<n-N_j)
			{
				r[i] = b[i] - aw[i]*x[i-N_j] - as[i]*x[i-1] - ap[i]*x[i] - an[i]*x[i+1] - ae[i]*x[i+N_j];
			}
			
			///from the end of east coefficient to the till end of north coefficient 444444444444444444444444444444444
			else if(i>=n-N_j && i<n-1)
			{
				r[i] = b[i] - aw[i]*x[i-N_j] - as[i]*x[i-1] - ap[i]*x[i] -as[i]*x[i+1];
			}
			
			//last row of A matrix
			else if(i==n-1)
			{
				r[i] = b[i] - aw[i]*x[i-N_j] - as[i]*x[i-1] - ap[i]*x[i];
			}
			
		}
		
		for(int i=0;i<n;i++)
		{
			//printf("r[%d]=%lf\n",i,r[i]);
		}
		
}
