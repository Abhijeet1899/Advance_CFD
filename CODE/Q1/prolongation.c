#include<stdio.h>

void prolongation(int m, double dh[2*m + 1], double d2h[m +1]  ) // a is at fine grid , b is at coarse grid
{

	for(int j=0;j<=m;j++)
	{
		dh[2*j]=d2h[j];// direct injection
		dh[2*j+1] = 0.5*(d2h[j] + d2h[j+1]); // interpolation
		//printf("sbdfhs\n");
		
	}
	
}
