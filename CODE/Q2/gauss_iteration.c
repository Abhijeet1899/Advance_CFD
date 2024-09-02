#include<stdio.h>

double gauss_iteration(int n, int neu, double x[n+1], double b[n+1]) // solving Ax=b equation using gauss siedel iterations for neu iterations
{
	double h=1.0/n;
	int k=0;
	double sum; 
	do{
		
		for(int i=1;i<n;i++)
		{
			x[i] = (x[i+1] + x[i-1] + b[i]*(h*h))/(2.0+h*h); 
		}
		
		k=k+1;
	}while(k<neu);
	
	//printf("gauss\n");
	
		
}
