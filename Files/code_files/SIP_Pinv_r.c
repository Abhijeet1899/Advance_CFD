#include<stdio.h>

double SIP_Pinv_r(int n,int N_j, double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double Pinvxr[n])
{
	double lp[n],lw[n],ls[n];
	double un[n],ue[n];
	double alpha = 0.5;
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
		if(i==0)
		{
			lp[i]=ap[i];
			un[i]=an[i]/lp[i];
			ue[i]=ae[i]/lp[i];
			ls[i]=as[i];
			lw[i]=aw[i];
		}
		
		else if(i>0 && i<n-N_j)
		{
			lw[i]=aw[i];
			ls[i]=as[i];
			lp[i]=ap[i] + ls[i]*(alpha*ue[i-1] + un[i-1]);
			un[i]=an[i]/lp[i];
			ue[i]=ae[i]/lp[i];
		}
		else if(i>=n-N_j && i<n)
		{
			lw[i]=aw[i]/(1.0+alpha*un[i-N_j]);
			ls[i]=as[i]/(1.0+alpha*ue[i-1]);
			lp[i]=ap[i]+ls[i]*(alpha*ue[i-1]+un[i-1])+lw[i]*(alpha*un[i-N_j]-ue[i-N_j]);
			un[i]=(an[i]-alpha*lw[i]*un[i-N_j])/lp[i];
			ue[i]=(ae[i]-alpha*ls[i]*ue[i-1])/lp[i];
		}
		
	}
	
	for(int j=0;j<n;j++)
	{
		//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",j+1,lw[j],ls[j],lp[j],un[j],ue[j]);
	}
	
	// solving LUx=r as its equal to Pinvxr=x; and P=LU
	
	double x[n],y[n];
	
	
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
