#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "head.h"

struct component {
	
		int m;
		double *v;
		double *f;

		
	};
	
	
	


int main()
{
	int n=512;
	int k=1;
	int neu=2;
	double c;
	double rh[n+1];
	c=pow(M_PI*k,2)+1.0;
	int l;
	l=n;
	int level=log(n)/log(2) -1;// no. of levels
	printf("level=%d\n",level);
	
	struct component F[level+1];
	
	
	
	for(int i=0;i<=level;i++)
	{
		F[i].m=l;
		printf("level%d n=%d\n",i+1,F[i].m);
		l=l/2;
	}
	
	
	
	
	for(int i=0;i<=level;i++)
    	{
    		F[i].v = (double *)malloc((F[i].m+1) * sizeof(double));
    		F[i].f = (double *)malloc((F[i].m+1)* sizeof(double));
    
 	   	
 	   	 if(F[i].v == NULL || F[i].f == NULL) {
            printf("Memory allocation failed1!\n");
            return 1;
    		}
 	}
 	
 	
 	  
    // initial guess for v and f
    for(int i=0;i<=level;i++)
    {
    for (int j = 0; j <= F[i].m; j++)
        {
            F[i].v[j]=0.0;
            F[i].f[j]=0.0;
            
        }
        
    }
    // initiallizing f
 	for (int j = 0; j <=n; j++)
        {
            F[0].f[j]  = c*sin(k*M_PI*j/n);
       
            
        }    
    printf("\n"); 
    
    // taking f to coaser grid
    for(int j=1;j<=level;j++)
    {
    	ristriction(F[j].m,F[j-1].f,F[j].f);
    
    }
   
    // applying v cycle level wise as discuss in FMG method
    for(int j=level;j>0;j--)
    {
    	gauss_iteration(F[j].m,neu,F[j].v,F[j].f);
    	
    	prolongation(F[j].m,F[j-1].v,F[j].v);
    	
    	v_cycle (F[j-1].m,  k, 1.0,1, F[j-1].v, F[j-1].f);
    	//v_cycle (F[j-1].m,  k, 1.0,1, F[j-1].v, F[j-1].f);
    	//v_cycle (F[j-1].m,  k, 1.0,1, F[j-1].v, F[j-1].f);
    	
    	
    }
    printf("FMG completed\n\n\n");
    // applying v cycle after FMG 
    v_cycle (F[0].m,  k, 1.0,2, F[0].v, F[0].f);
    residual(F[0].m,rh,F[0].v,F[0].f);
    
    for (int j = 0; j<=n; j++)
        { 
           //printf("vh[%d]=%lf\n",j, F[0].v[j]);
           //printf("u[%d]=%lf\n",j,sin(k*M_PI*j/n));
           //printf("rh[%d]=%lf\n",j,rh[j]);
           
        }
        
    //printf("\n");
    
    //v_cycle(n,k,1.0,F[0].v,F[0].f);
   /* for(int i = 0; i <= level; i++) {
        free(F[i].v);
        free(F[i].f);
        
    }*/
     
     
	
	return 0;	
 }
 
 
 
 

