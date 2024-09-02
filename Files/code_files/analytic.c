#include<stdio.h>
#include<math.h>


int main()
{
	int n=4;
	double sum;
	double t[n][n];
	double x,y;
	x=0.0;y=0.0;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			sum=0.0;
			for(int k=1;k<50;k++)
			{
				sum=sum+((pow(-1,1+k)+1)/k)*sin(k*M_PI*x)*(sinh(k*M_PI*y)/sinh(k*M_PI));
			}
			t[i][j]=(2.0/M_PI)*sum;
			y=y+j/128;
			
		}
		x=x+i/128;
	}
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			printf("T[%d][%d]=%lf  ",j,i,t[j][i]);
		}
		printf("\n");
	}
	
}
