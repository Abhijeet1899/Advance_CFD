#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "head.h"




int main()
{
	int N,N_i,N_j;
	N_i=128;N_j=128;
	N=N_i*N_j;
	double ap[N],aw[N],ae[N],an[N],as[N];// A coefficent matrix divided in 5 column matrix
	double b[N];// right hand side vector
	double x[N];// temprature vector
	double r[N],r_pre[N];// residual vector
	double d[N]; // direction vector
	double Ad[N];///vector to store multiplication of A matrix with d (direction) vector 
	double alpha, beta, g,q,epsilon=0.000001; //here g is the multiplication of r tvector transpose and r
	int i,j,k,l;
	// Ax=b;
	i=2;
	j=1;
	
	FILE *file;
	file=fopen("result.txt","w");
	//initializing A sparse matrix using 5 colum matrix of size N (ap,aw,as,an,ae)
	for(int n=0;n<N;n++)
	{
		// bottom left corner i=0
		if(n==0){
			ap[n]=6.0;
			ae[n]=-1.0;an[n]=-1.0;
			as[n]=0.0;aw[n]=0.0;
			b[n]=0.0;
		}	
		
		//left boundary T=0 at boundary
		else if(n>0 && n<N_j-1){
			ap[n]=5.0;
			ae[n]=-1.0;an[n]=-1.0;as[n]=-1.0;
			aw[n]=0.0;
			b[n]=0.0;
		}
		
		//left top corner
		else if(n==N_j-1){
			ap[n]=6.0;
			ae[n]=-1.0;as[n]=-1.0;
			an[n]=0.0;aw[n]=0.0;
			b[n]=2.0;
		}
		
		//top boundary T=1
		else if(n==i*N_j-1 && i<N_i)
		{
			ap[n]=5.0;
			ae[n]=-1.0;aw[n]=-1.0;as[n]=-1.0;
			an[n]=0.0;
			b[n]=2.0;
			//printf("i=%d\t",i);
			i=i+1;
			//printf("n=%d\n",n);
		}
		
		//top right corner
		else if(n==N-1){
			ap[n]=6.0;
			aw[n]=-1.0;as[n]=-1.0;
			an[n]=0.0;ae[n]=0.0;
			b[n]=2.0;
			//printf("ntr=%d\n",n);
		}
		
		//right boundary T=0
		else if(n>((N_i-1)*N_j) && n<=N_i*N_j-2){
			ap[n]=5.0;
			aw[n]=-1.0;an[n]=-1.0;as[n]=-1.0;
			ae[n]=0.0;
			b[n]=0.0;
			//printf("nrb=%d\n",n);
		}
		
		//bottom right corner
		else if(n==((N_i-1)*N_j)){
			ap[n]=6.0;
			aw[n]=-1.0;an[n]=-1.0;
			ae[n]=0.0;as[n]=0.0;
			b[n]=0.0;
			printf("\n n=%d\n",n);
		}
		
		//bottom boundary T=0
		else if(n==j*N_j && j<N_i){
			ap[n]=5.0;
			ae[n]=-1.0;aw[n]=-1.0;an[n]=-1.0;
			as[n]=0.0;
			b[n]=0.0;
			//printf("j=%d\t",j);
			j=j+1;
			//printf("nbb=%d\n",n);
		}
		
		// interior points
		else{
			ap[n]=4.0;
			ae[n]=-1.0;an[n]=-1.0;as[n]=-1.0;aw[n]=-1.0;
			b[n]=0.0;
		}
			
	}
	
	
	/// initiallizing the temprature (here Temp is taken x)
	for(int n=0;n<N;n++)
	{
		x[n]=0.0;
	}
	
	residual(N,N_j,ap,aw,as,an,ae,r,x,b);// finding rsidual
	
	for(int n=0;n<N;n++)
	{
		//printf("aw[%d]=%lf  as[%d]=%lf  ap[%d]=%lf  an[%d]=%lf  ae[%d]=%lf x[%d]=%lf b[%d]=%lf r[%d]=%lf\n",
		//n,aw[n],n,as[n],n,ap[n],n,an[n],n,ae[n],n,x[n],n,b[n],n,r[n]);
	}
	
	// for first search direction is the residual vector itself
	for(int n=0;n<N;n++)
	{
		d[n]=r[n];
	}
	
	// g is the sacalar to store vector product
	g=vectorDotproduct(N,r,r);
	
	//to get values for vector after multiplication of A and direction vector d
	matrix_vector_product(N,N_j,ap,aw,as,an,ae,d,Ad);
	
	// calculating alpha
	alpha=g/vectorDotproduct(N,d,Ad);
	
	//printf("alpha=%lf\n",alpha);
	k=0;// iteration dialier
	l=1;// for calculating residual using r=b-Ax
	
	// iterative loop
	do
	{
		for(int n=0;n<N;n++)
		{
			x[n]=x[n]+alpha*d[n]; // updating temprature (using directional search
		}
		
		if(k==l*50){
			residual(N,N_j,ap,aw,as,an,ae,r,x,b); // finding residual
			l=l+1;
			//printf("l=%d\n",l);
		}
		else
		{
			for(int n=0;n<N;n++)
			{
				r[n]=r[n]-alpha*Ad[n]; // finding residual
			}
			//printf("////////\n");
		}
		
		beta=-vectorDotproduct(N,r,r)/g; // calculating beta
		
		for(int n=0;n<N;n++)
			{
				d[n]=r[n]-beta*d[n]; // updating search direction
			}
		
		// calculating values that will be used in next iteration
		g=vectorDotproduct(N,r,r);
		matrix_vector_product(N,N_j,ap,aw,as,an,ae,d,Ad);
		alpha=g/vectorDotproduct(N,d,Ad);
		//printf("alpha=%lf\n",alpha);
		k=k+1;
		
		
		q=two_norm(N,r);/// calculating 2nd norm of residual for checking termination condition
		//printf("q=%lf\n",q);
		fprintf(file,"%d\t%lf\n",k,q);
		//printf("iterations=%d\n",k);
			
	}while(k<1000 && q>epsilon);
		
		printf("iterations=%d\n",k);
		printf("q=%lf\n",q);
		for(int n=0;n<N;n++)
			{
				//printf("residual=%lf\n",r[n]);
			}
		
	
	
	
	
}
