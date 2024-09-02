#include<stdio.h>

void ristriction(int m, double dh[2*m+1], double d2h[m+1]  )// a is at fine grid , b is at coarse grid
{
	d2h[0]=dh[0];
	d2h[m]=dh[2*m];
	// direct injection
	for(int j=1;j<m;j++)
	{
		//d2h[j]=dh[2*j]; // direct injection
		
		d2h[j]=(dh[2*j-1] + 2.0*dh[2*j] + dh[2*j+1])/4.0;
		
	}
	
}
