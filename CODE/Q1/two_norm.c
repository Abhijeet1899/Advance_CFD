#include<stdio.h>
#include<math.h>

double two_norm(int n, double r[n+1])
{
	double d;
	
		d=0.0;
		for(int j=0;j<=n;j++)
		{
			
				d=d+pow(r[j],2);
		}
	
	return sqrt(d/(n+1));
	 
}
