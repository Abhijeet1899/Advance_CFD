#include<stdio.h>
#include<math.h>



double vectorDotproduct(int n,double vec1[n],double vec2[n]) //// vector dot product
{
		double result,sum=0.0;
		for(int j=0;j<n;j++)
		{
			sum=sum+vec1[j]*vec2[j];
		}	
		result=sum;
		//	printf("result=%lf",result);
			return result;
		
}   
