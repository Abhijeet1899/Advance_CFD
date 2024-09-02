#include<stdio.h>

double ILU_Pinv_r(int n,int N_j, double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double Pinvxr[n])
{
	double lp[n],lw[n],ls[n];
	double un[n],ue[n];
	double x[n],y[n];
	
	for(int i=0;i<n;i++)
	{
		lp[i]=0.0;
		lw[i]=0.0;
		lp[i]=0.0;
		un[i]=0.0;
		ue[i]=0.0;
	}
	// calculating L and U matrix (column matrices )
	for(int i=0;i<n;i++)
	{
		lw[i]=aw[i];
		ls[i]=as[i];
		
		
		if(i==0){
			lp[i]=ap[i];
		}
		else if (i>=1 && i<=N_j-1){
			lp[i]=ap[i]-ls[i]*un[i-1];
		}
		else{
			lp[i]=ap[i]-ls[i]*un[i-1]-lw[i]*ue[i-N_j+1];
		}
		
		un[i]=an[i]/lp[i];
		ue[i]=ae[i]/lp[i];
		
	}
	
	for(int j=0;j<n;j++)
	{
		//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",j+1,lw[j],ls[j],lp[j],un[j],ue[j]);
	}
	
	// solving LUx=r as its equal to Pinvxr=x; and P=LU
	
	
	
	
	/// solving Ly=r and getting y , here y is Ux
	//forward substitution
	for(int i=0;i<n;i++)
	{
		if(i==0){
			y[i]=r[i]/lp[i];
			//printf("y[%d]=%lf\n",i,y[i]);
		}
		/*else if(i>0 && i<=N_j-1){
			y[i]=(r[i]-ls[i]*y[i-1])/lp[i];
			//printf("y[%d]=%lf\n",i,y[i]);
		}*/
		else{
			y[i]=(r[i]-ls[i]*y[i-1]-lw[i]*y[i-N_j])/lp[i];
			//printf("y[%d]=%lf\n",i,y[i]);
		}
	}
	
	
	//solving Ux=y , getting x, which is multiplication of P inverse and r , gives Pinvxr
	// backward substitution
	
	for(int i=n-1;i>=0;i--)
	{
		if(i==n-1){
			x[i]=y[i];
		}
		else if(i<n && i>=n-N_j){
			x[i]=y[i]-un[i]*x[i+1];
		}
		else{
			x[i]=y[i]-un[i]*x[i+1]-ue[i]*x[i+N_j];
		}
	}
	
	for(int i=0;i<n;i++)
	{
		Pinvxr[i]=x[i];
		//printf("Pinvxr[%d]=%lf\n",i,Pinvxr[i]);
	}
	
	return 0;
	
	
}
