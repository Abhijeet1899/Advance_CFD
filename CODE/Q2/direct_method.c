#include<stdio.h>

double gauss_elemination(int n, double A[][n],double c[n], double x[n])
{
	double agu[n][n+1];
	double a[n][n];
	double b[n];
	double sum, constant,num;
	//gauss elimination
	
	 for(int i=0;i<n;i++)
    	{
        	for(int j=0;j<n;j++)
        	{
        	   agu[i][j]=A[i][j];
        	}
       	
    	}
    	for(int j=0;j<n;j++)
        	{
        	   agu[j][n]=c[j];
        	}
    	
	//zero pivoting problem
	for(int j=0;j<n+1;j++)
	{
		for(int i=j;i<n;i++)
		{
			if(agu[j][j]==0 && agu[i][j]!=0)               //swap the rows
			{
				for(int k=0;k<n+1;k++)
				{
					num=agu[j][k];
					agu[j][k]=agu[i][k];
					agu[i][k]=num;
				}
			}
		}
	}
		for(int j=0;j<n+1;j++)
		{
			for(int i=0;i<n;i++)
			if(i>j)
			{
				constant=agu[i][j]/agu[j][j];
				for(int k=0;k<n+1;k++)
				{
					agu[i][k]=agu[i][k]-constant*agu[j][k];
				}
			}
			
		
		}
		
	// extracting echlon for coefficient matrix and column matrix
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			
			a[i][j]=agu[i][j];
		}
		
			
		
			b[i]=agu[i][n];
		
	}
	
	// back substitution for getting the solution of the linear equations
	x[n-1]=b[n-1]/a[n-1][n-1];
	
	for(int i=n-2;i>=0;i--)
	{
		sum=0;
		for(int j=i+1;j<n;j++)
		{
			sum=sum+a[i][j]*x[j];
		}
		x[i]=(b[i]-sum)/a[i][i];
		
	}
	
	//printing the solutions
	for(int i=0;i<n;i++)
	{
		printf("x%d=%lf   ",i,x[i]);
	}
	
	
}
