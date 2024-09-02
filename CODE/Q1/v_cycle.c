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
	

double v_cycle (int n, int k )
{
	
	int l;
	int level=log(n)/log(2) -1;// no. of levels
	//printf("level=%d\n",level);
	int neu=2;/// no. of gauss iterations in each level
	l=n;
	
	double temp;
	double c;
	c=pow(M_PI*k,2)+1.0;
	
	struct element v[level+1];
	
	FILE *file;
	file=fopen("plot.txt","w");
	
	//inserting no. of grid point in each level
	for(int i=0;i<=level;i++)
	{
		v[i].m=l;
		//printf("level (%d) n=%d\n",i+1,v[i].m);
		l=l/2;
	}
	
	//printf("level n=%d\n",v[level].m);
	
	//assigning memory for the variables at each level
	for(int i=0;i<=level;i++)
    	{
    		v[i].v = (double *)malloc((v[i].m+1) * sizeof(double));
    		v[i].f = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	v[i].r = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	v[i].e = (double *)malloc((v[i].m+1)* sizeof(double));
 	   	
 	   	 if(v[i].r == NULL || v[i].e == NULL || v[i].v == NULL || v[i].f == NULL) {
            printf("Memory allocation failed!\n");
          
    		}
    
 	}
	
	
    
    
    // initial guess and assigning all others to zero initially
    for(int i=0;i<=level;i++)
    {
    for (int j = 0; j <= v[i].m; j++)
        {
            v[i].v[j]=0.0;
            v[i].e[j]=0.0;
            v[i].r[j]=0.0;
            v[i].f[j]=0.0;
        }
        
    }
    
   // assigning values to f in Av=f equation at fine grid 
   for (int j = 0; j <=n; j++)
        {
            v[0].f[j]  = c*sin(k*M_PI*j/n);
       
            
        }    
    
     
 	
 	int s=0;
 do{	
 
 for(int i=1;i<=level;i++)
    {
    for (int j = 0; j <= v[i].m; j++)
        {
            v[i].v[j]=0.0;
            v[i].e[j]=0.0;
            v[i].r[j]=0.0;
            v[i].f[j]=0.0;
        }
        
    }
 
   /// v cycle down operation 
    for(int i=0;i<level;i++)
    {
    	gauss_iteration(v[i].m,neu,v[i].v,v[i].f);
    	
    	residual(v[i].m,v[i].r,v[i].v,v[i].f);
    	
    	ristriction(v[i+1].m,v[i].r,v[i+1].f);
    	
    	 
    }
    
    gauss_iteration(v[level].m,neu,v[level].v,v[level].f);////direct method solution at the last coarser grid
    
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
    
    // finding rsidual at the end of v cycle
    residual(v[0].m,v[0].r,v[0].v,v[0].f);
    
    // 2nd norm of the residual
    temp=two_norm(n,v[0].r);
    //printf("iteration=%d\t temp=%lf\n",s,temp);
    fprintf(file,"%d\t%lf\n",s+1,temp);
    s=s+1;
    }while(temp>0.000001 && s<1500);
    
    printf("iteration=%d\n",s);
    
    for (int j = 0; j<=n; j++)
        { 
          /* printf("vh[%d]=%lf\n",j, v[0].v[j]);
           printf("u[%d]=%lf\n",j,sin(k*M_PI*j/n));
           printf("rh[%d]=%lf\n",j,v[0].r[j]);
           */
        }
        
  
    
   
    
 }
