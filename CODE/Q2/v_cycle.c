#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "head.h"
struct element {
	
		int m;
		double *v;
		double *f;
		double *r;
		double *e;
		 
	};


double v_cycle (int n, int k, double sigma,int a, double vh[n+1], double fh[n+1])
{
	int l;
	int level=log(n)/log(2) -1;// no. of levels
	printf("V_cycle-level=%d\n",level);
	int neu=2;/// no. of gauss iterations in each level
	l=n;
	double temp;
	double c;
	c=pow(M_PI*k,2)+sigma;
	
	struct element v[level+1];
	

	for(int i=0;i<=level;i++)
	{
		v[i].m=l;
		printf("v_cycle-levels%d n=%d\n",i,v[i].m);
		l=l/2;
	}
	
	for(int i=0;i<=level;i++)
    	{
    		v[i].v = (double *)malloc((v[i].m+1) * sizeof(double));
    		v[i].f = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	v[i].r = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	v[i].e = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	
 	   	 if(v[i].r == NULL || v[i].e == NULL || v[i].v == NULL || v[i].f == NULL) {
            printf("Memory allocation failed!\n");
            return 1;
    		}
 	}
 	
 	 
	
    
    // initial guess for v and error
    for(int i=0;i<=level;i++)
    {
    for (int j = 0; j <=v[i].m; j++)
        {
            v[i].v[j]=0.0;
            v[i].f[j]=0.0;
            v[i].r[j]=0.0;
            v[i].e[j]=0.0;
            
        }
        
    }
    
    // input from FMG
   for (int j = 0; j <=n; j++)
        {
            //v[0].f[j]  = c*sin(k*M_PI*j/n);
            
           v[0].f[j] = fh[j];
            v[0].v[j] = vh[j];
       
            
        }    
    
     
 	
 	int s=0;
 	
 do{
 
 	 for(int i=1;i<=level;i++)
    {
    for (int j = 0; j <=v[i].m; j++)
        {
            v[i].v[j]=0.0;
            v[i].f[j]=0.0;
            v[i].r[j]=0.0;
            v[i].e[j]=0.0;
            
        }
        
    }	
   /// v cycle down operation 
    for(int i=0;i<level;i++)
    {
    	gauss_iteration(v[i].m,neu,v[i].v,v[i].f);
    	
    	residual(v[i].m,v[i].r,v[i].v,v[i].f);
    	
    	ristriction(v[i+1].m,v[i].r,v[i+1].f);
    	
    	 
    }
    
    gauss_iteration(v[level].m,neu,v[level].v,v[level].f);////
    
    ///v cycle up operation
    for(int i=level;i>0;i--)
    {
    	prolongation(v[i].m,v[i-1].e,v[i].v);
    	
    	 for(int j=0;j<=v[i-1].m;j++)
   	 {
    	
    		v[i-1].v[j] = v[i-1].v[j] + v[i-1].e[j];
   	 }
   	 	
   	  gauss_iteration(v[i-1].m,neu,v[i-1].v,v[i-1].f);
   	  	 
    }
    
    residual(v[0].m,v[0].r,v[0].v,v[0].f);
    
    temp=two_norm(n,v[0].r);
    printf("temp=%lf\n",temp);
    s=s+1;
    }while((a==1 && s<2) || (a==2 && temp>0.000001));
    
    printf("\nV-cycle-iteration=%d\n",s);
    for (int j = 0; j<=n; j++)
        { 
           //printf("vh[%d]=%lf\n",j, v[0].v[j]);
           //printf("u[%d]=%lf\n",j,sin(k*M_PI*j/n));
           //printf("rh[%d]=%lf\n",j,v[0].r[j]);
           
           fh[j]= v[0].f[j] ;
           vh[j]= v[0].v[j] ;
           
        }
        
   // printf("\n");
    
   // Free allocated memory
/*for(int i = 0; i <= level; i++) {
    free(v[i].v);
    free(v[i].f);
    free(v[i].r);
    free(v[i].e);
}*/
    return 0;
    
 }
 
 
 
