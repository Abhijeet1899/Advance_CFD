#include<stdio.h>

double residual(int n,double r[n+1], double x[n+1], double b[n+1])
{
	double h=1.0/n;
	
		for(int i=0;i<=n;i++)
		{
			if(i>0 && i<n)
			{
				r[i] = b[i] -( (-x[i+1] - x[i-1] + 2.0*x[i])/(h*h) + x[i]); 
			} 
			else{r[i]=0.0;}
		}
		
}
