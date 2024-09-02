#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// function for applying SIMPLE method
void SIMPLE (int n, double D, double h, double Re, double p_star[][n+1],double p_dash[][n+1], 
			double u_star[n+1][n+1], double v_star[n+1][n+1], double a[2][n+1][n+1], 
			double aw[2][n+1][n+1], double ae[2][n+1][n+1],double as[2][n+1][n+1], double an[2][n+1][n+1], 
			double fw[2][n+1][n+1], double fe[2][n+1][n+1],double fs[2][n+1][n+1], double fn[2][n+1][n+1] ,
			double sp[2][n+1][n+1], double su[2][n+1][n+1]) ;
			
// function for updating the coefficient of the momentum equations
void updating_coefficient(int n, double D, double u_star[n+1][n+1], double v_star[n+1][n+1], double a[2][n+1][n+1], 
				double aw[2][n+1][n+1], double ae[2][n+1][n+1],double as[2][n+1][n+1], double an[2][n+1][n+1] ,
				double fw[2][n+1][n+1], double fe[2][n+1][n+1],double fs[2][n+1][n+1], double fn[2][n+1][n+1], 
				double sp[2][n+1][n+1], double su[2][n+1][n+1]); 

// function for solving momemtum equation						
void momentum_eq(int n, double D, double h, double p_star[][n+1], double u_star[n+1][n+1], double v_star[n+1][n+1], double a[2][n+1][n+1], 
		double aw[2][n+1][n+1], double ae[2][n+1][n+1],double as[2][n+1][n+1], double an[2][n+1][n+1],
		 double fw[2][n+1][n+1], double fe[2][n+1][n+1],
				double fs[2][n+1][n+1], double fn[2][n+1][n+1], 
				double sp[2][n+1][n+1], double su[2][n+1][n+1]) ;

//function for solving pressure correction equation						
void p_dash_equation(int n, double h, double u_star[n+1][n+1], double v_star[n+1][n+1], double p_dash[][n+1], double a[2][n+1][n+1]);

// function for applying hybrid scheme (max operator)
double max(double i, double j);


// function for printing values
void print_onscreen(int n, double g[n+1][n+1]);

double rms_termination_condition(int n, double u[n+1][n+1], double uold[n+1][n+1]);// 2nd norm 
double infinite_norm(int n, double u_star[n+1][n+1], double uold[n+1][n+1]);// infinity norm


int main()
{
	int n=128;
	double l;
	double h=1.0/n;
	double Re=100.0;
	double D=(1.0/(h*Re)) ;
	double u[n+1][n+1], u_star[n+1][n+1], v[n+1][n+1], v_star[n+1][n+1], p[n+1][n+1], p_star[n+1][n+1], p_dash[n+1][n+1], b[n+1][n+1];
	double a[2][n+1][n+1];
	double aw[2][n+1][n+1],ae[2][n+1][n+1], as[2][n+1][n+1],an[2][n+1][n+1];
	double  fw[2][n+1][n+1], fe[2][n+1][n+1], fs[2][n+1][n+1], fn[2][n+1][n+1];	
	double sp[2][n+1][n+1], su[2][n+1][n+1];
	double alpha_p, alpha_u, alpha_v;
	double w[n+1][n+1];
	double psi[n+1][n+1];
	// initial values
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			
				u_star[i][j]=0.0;
				v_star[i][j]=0.0;
				p_star[i][j]=0.0;
				su[0][i][j]=0.0;
				su[1][i][j]=0.0;
				sp[0][i][j]=0.0;
				sp[1][i][j]=0.0;
				
				w[i][j] = 0.0;
				psi[i][j] = 0.0;
		}
	}	
	
	// applying SIMPLE method
	SIMPLE (n,D,h,Re,p_star,p_dash,u_star,v_star,a,aw,ae,as,an,fw,fe,fs,fn,sp,su );	
	
	
	
	// solving vorticity using known velocities
	for (int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			w[i][j] = (v_star[i][j] - v_star[i][j-1] - u_star[i][j] + u_star[i-1][j])/h;
		}
	}
	
	
	

	
	// printing results in file
	double x=0.0,y=0.0;
	FILE *fp;
	fp=fopen("output_100.dat","w");
	fprintf(fp,"ZONE I=%d, J=%d\n",n+1,n+1);
	for(int i=0;i<=n;i++)
	{
		y=i*h;
		for(int j=0;j<=n;j++)
		{
			x=j*h;
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n",x,y,u_star[i][j],v_star[i][j],w[i][j],psi[i][j]);
		}
		
	}
	fclose(fp); 

	
	FILE *f;
	f=fopen("u_mid_100.dat","w");
	//fprintf(f,"ZONE I=%d, J=%d\n",m,n);
	
		for(int j=0;j<=n+1;j++)
		{
			if(j==0){x=0.0;fprintf(f,"%lf\t%lf\n",u_star[j][64],x);}
			else if(j==1){x=0.5*h*j;fprintf(f,"%lf\t%lf\n",u_star[j][64],x);}
			else if(j>1 && j<n+1 ){
			x=h*j-0.5*h;
			fprintf(f,"%lf\t%lf\n",u_star[j][64],x);	
			}
			else
			{
				x=h*(j-1);
			fprintf(f,"%lf\t%lf\n",1.0,x);
			}
			
		
	}
	fclose(f);
	
	FILE *f1;
	f1=fopen("v_mid_100.dat","w");
	//fprintf(f,"ZONE I=%d, J=%d\n",m,n);
	
		for(int j=0;j<=n;j++)
		{
			x=h*j;
			
			
			
				fprintf(f,"%lf\t%lf\n",v_star[64][j],x);	
			
		
	}
	fclose(f1);
	
	
	
}


// defining SIMPLE function
void SIMPLE (int n, double D, double h, double Re, double p_star[][n+1],double p_dash[][n+1], 
			double u_star[n+1][n+1], double v_star[n+1][n+1], double a[2][n+1][n+1], 
			double aw[2][n+1][n+1], double ae[2][n+1][n+1],double as[2][n+1][n+1], double an[2][n+1][n+1], 
			double fw[2][n+1][n+1], double fe[2][n+1][n+1],double fs[2][n+1][n+1], double fn[2][n+1][n+1] ,
			double sp[2][n+1][n+1], double su[2][n+1][n+1]) 
{
	double e1,e2,e3,epsilon;
	e1=0.0;e2=0.0;e3=0.0;
	epsilon=0.00001;
	double u_old[n+1][n+1];
	double v_old[n+1][n+1];
	double pstar_old[n+1][n+1];
	
	
	/// Initial Guess 
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			
				u_star[i][j]=0.0;
				v_star[i][j]=0.0;
				p_star[i][j]=0.0;
				su[0][i][j]=0.0;
				su[1][i][j]=0.0;
				sp[0][i][j]=0.0;
				sp[1][i][j]=0.0;
				u_old[i][j]=0.0;
				v_old[i][j]=0.0;
				pstar_old[i][j]=0.0;
				
		}
	}
	
	
	int k=0;//to keep track of no. of SIMPLE iterations
	
	
	do{
		
	
		
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			
				
				p_dash[i][j]=0.0;
				
		}
	}
	
	
	
	
	
		
		
		momentum_eq(n,D,h,p_star,u_star,v_star,a,aw,ae,as,an,fw,fe,fs,fn,sp,su );  /// solving momentum equation
		
		
		
		p_dash_equation(n, h,u_star, v_star,p_dash,a); // solving correction pressure equation (by using velocity in terms of pressure correction and using continuity equation
		
		
		
		
	
	
	//// Updating pressure term
		for(int i=1;i<=n;i++)
		{
			for(int j=1;j<=n;j++)
			{
			
				p_star[i][j]=p_star[i][j]+p_dash[i][j];
				
				
				
			}
		}
		
		///Updating u velocity term
		for(int i=1;i<=n;i++)
		{
			for(int j=1;j<=n-1;j++)
			{
			
				
				u_star[i][j]=u_star[i][j] + (h/a[0][i][j])*(p_dash[i-1][j]-p_dash[i][j]);
				
			}
		}
		
		/// updating v velocity term
		for(int i=1;i<=n-1;i++)
		{
			for(int j=1;j<=n;j++)
			{
			
				
				
				v_star[i][j]=v_star[i][j] + (h/a[1][i][j])*(p_dash[i][j-1]-p_dash[i][j]);
				
			}
		}
		
		
		/// updating u velocity with under relaxation factor 0.5
		for(int i=1;i<=n;i++)
		{
			for(int j=1;j<=n-1;j++)
			{
			
				
				u_star[i][j]=0.5*u_star[i][j] + 0.5*u_old[i][j];
				
			}
		}
		
		/// updating v velocity with under relaxation factor 0.5
		for(int i=1;i<=n-1;i++)
		{
			for(int j=1;j<=n;j++)
			{
			
				
				
				v_star[i][j]=0.5*v_star[i][j] + 0.5*v_old[i][j];
				
			}
		}
		
		// updating P-star using under relaxation factor 0.5
		for(int i=1;i<=n;i++)
		{
			for(int j=1;j<=n;j++)
			{
			
				
				p_star[i][j]=0.5*p_star[i][j] + 0.5*pstar_old[i][j];
				
			}
		}
		
		
		
		e1=rms_termination_condition(n, u_star,u_old); // 2nd norm of u velocity
		e2=rms_termination_condition(n, v_star,v_old); // 2nd norm of v velocity
		e3=rms_termination_condition(n, p_star,pstar_old); // 2nd norm of pressure
		
		
		// storing previous iteration velocities
		for(int i=0;i<=n;i++)
		{
			for(int j=0;j<=n;j++)
			{
				
					u_old[i][j] = u_star[i][j];
					v_old[i][j] = v_star[i][j];
					pstar_old[i][j] = p_star[i][j];
					
					
			}
		}	
		k=k+1;
		printf("k=%d\n",k);
		printf("e_u=%lf\t e_v=%lf\t e3=%lf\n",e1,e2,e3);
	}while (e1>epsilon || e2>epsilon );
	
	
	
	
	
	
	
}


/// function to evaluate coefficient of the momentum equation
void updating_coefficient(int n, double D, double u_star[n+1][n+1], double v_star[n+1][n+1], 
				double a[2][n+1][n+1], 
				double aw[2][n+1][n+1], double ae[2][n+1][n+1],
				double as[2][n+1][n+1], double an[2][n+1][n+1] ,
				double fw[2][n+1][n+1], double fe[2][n+1][n+1],
				double fs[2][n+1][n+1], double fn[2][n+1][n+1], 
				double sp[2][n+1][n+1], double su[2][n+1][n+1]) 
{
	/// boundary coeefficient updation for u momentum equation
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n-1;j++)
		{
		
		
		if (i==1 ) /// bottom boundary condition coeefficients with corners 
		{
			fw[0][i][j] = 0.5*(u_star[i][j]+u_star[i][j-1]);
			fe[0][i][j] = 0.5*(u_star[i][j+1]+u_star[i][j]);
			fs[0][i][j] = 0.0 ;// 0.5*(v_star[i-1][j+1]+v_star[i-1][j]);/// i -- >i+1, i-1 --> i
			fn[0][i][j] = 0.5*(v_star[i][j+1]+v_star[i][j]);		
			
			
			aw[0][i][j] = max(fw[0][i][j], D + 0.5*fw[0][i][j]);
			ae[0][i][j] =  max(-fe[0][i][j],D - 0.5*fe[0][i][j]);
			as[0][i][j] =  0.5*fs[0][i][j]; // Ds is zero
			an[0][i][j] =  max(-fn[0][i][j],D - 0.5*fn[0][i][j]);		
			sp[0][i][j] = -2.0*D;// variable source term from difussive term
			
			
			
			
			a[0][i][j] = aw[0][i][j] + ae[0][i][j] + as[0][i][j] + an[0][i][j] 
			+ fe[0][i][j] - fw[0][i][j] + fn[0][i][j] - fs[0][i][j] - sp[0][i][j];
			
		}
		
		
		
		else if(i==n ) /// top boundary condition with corners
		{
			fw[0][i][j] = 0.5*(u_star[i][j]+u_star[i][j-1]);
			fe[0][i][j] = 0.5*(u_star[i][j+1]+u_star[i][j]);
			fs[0][i][j] = 0.5*(v_star[i-1][j+1]+v_star[i-1][j]);
			fn[0][i][j] = 0.0;//0.5*(v_star[i][j+1]+v_star[i][j]);	// fn is zero	
			
			
			aw[0][i][j] = max(fw[0][i][j],D + 0.5*fw[0][i][j]);
			ae[0][i][j] =  max(-fe[0][i][j],D - 0.5*fe[0][i][j]);
			as[0][i][j] =  max(fs[0][i][j],D + 0.5*fs[0][i][j]); 
			an[0][i][j] =  0.5*fn[0][i][j];	// Dn is zero	
			sp[0][i][j] = -2.0*D; // variable source term from difussive term
			su[0][i][j]=2.0*D;  /// upper plate is moving we got source term 
			
			
			
			a[0][i][j] = aw[0][i][j] + ae[0][i][j] + as[0][i][j] + an[0][i][j] 
			+ fe[0][i][j] - fw[0][i][j] + fn[0][i][j] - fs[0][i][j] - sp[0][i][j];
			
		} 
		
		
		else if( j==1 && (i>1 && i<n)) // left boundary without corners
		{
			fw[0][i][j] = 0.5*(u_star[i][j]+u_star[i][j-1]);
			fe[0][i][j] = 0.5*(u_star[i][j+1]+u_star[i][j]);
			fs[0][i][j] = 0.5*(v_star[i-1][j+1]+v_star[i-1][j]);
			fn[0][i][j] = 0.5*(v_star[i][j+1]+v_star[i][j]);		
			
			
			aw[0][i][j] = max(fw[0][i][j],D + 0.5*fw[0][i][j]);
			ae[0][i][j] =  max(-fe[0][i][j],D - 0.5*fe[0][i][j]);
			as[0][i][j] = max(fs[0][i][j],D + 0.5*fs[0][i][j]);
			an[0][i][j] =  max(-fn[0][i][j],D - 0.5*fn[0][i][j]);		
			
			
			
			
			
			
			a[0][i][j] = aw[0][i][j] + ae[0][i][j] + as[0][i][j] + an[0][i][j] 
			+ fe[0][i][j] - fw[0][i][j] + fn[0][i][j] - fs[0][i][j] - sp[0][i][j];
			
		}
		else if( j==n-1 && (i>1 && i<n)) // right boundary
		{
			fw[0][i][j] = 0.5*(u_star[i][j]+u_star[i][j-1]);
			fe[0][i][j] = 0.5*(u_star[i][j+1]+u_star[i][j]);
			fs[0][i][j] = 0.5*(v_star[i-1][j+1]+v_star[i-1][j]);
			fn[0][i][j] = 0.5*(v_star[i][j+1]+v_star[i][j]);		
			
			
			aw[0][i][j] = max(fw[0][i][j],D + 0.5*fw[0][i][j]);
			ae[0][i][j] = max(-fe[0][i][j], D - 0.5*fe[0][i][j]);
			as[0][i][j] = max(fs[0][i][j],D + 0.5*fs[0][i][j]);
			an[0][i][j] =  max(-fn[0][i][j],D - 0.5*fn[0][i][j]);		
			
			
			
			
			
			
			a[0][i][j] = aw[0][i][j] + ae[0][i][j] + as[0][i][j] + an[0][i][j] 
			+ fe[0][i][j] - fw[0][i][j] + fn[0][i][j] - fs[0][i][j] - sp[0][i][j];
			
		}
	}	
	}
	
	//u-momentum coefficient for interior points
	for(int i=2;i<=n-1;i++)
	{
		for(int j=2;j<=n-2;j++)
		{
		
		
			fw[0][i][j] = 0.5*(u_star[i][j]+u_star[i][j-1]);
			fe[0][i][j] = 0.5*(u_star[i][j+1]+u_star[i][j]);
			fs[0][i][j] = 0.5*(v_star[i-1][j+1]+v_star[i-1][j]);
			fn[0][i][j] = 0.5*(v_star[i][j+1]+v_star[i][j]);		
			
			
			aw[0][i][j] = max(fw[0][i][j],D + 0.5*fw[0][i][j]);
			ae[0][i][j] = max(-fe[0][i][j], D - 0.5*fe[0][i][j]);
			as[0][i][j] = max(fs[0][i][j],D + 0.5*fs[0][i][j]);
			an[0][i][j] =  max(-fn[0][i][j],D - 0.5*fn[0][i][j]);		
			
			
			
			
			
			
			a[0][i][j] = aw[0][i][j] + ae[0][i][j] + as[0][i][j] + an[0][i][j] 
			+ fe[0][i][j] - fw[0][i][j] + fn[0][i][j] - fs[0][i][j] - sp[0][i][j];
			
			
				
		}
	}
	
	
	
	
	// y- momentum equations boundary
		for(int i=1;i<=n-1;i++)
	{
		for(int j=1;j<=n;j++)
		{
		
			if (j==1) //left boundary
			{
			
			fw[1][i][j] = 0.0;// 0.5*(u_star[i+1][j-1]+u_star[i][j-1]);
			fe[1][i][j] = 0.5*(u_star[i+1][j]+u_star[i][j]);
			fs[1][i][j] = 0.5*(v_star[i][j]+v_star[i-1][j]);
			fn[1][i][j] = 0.5*(v_star[i+1][j]+v_star[i][j]);
			
		
			
			aw[1][i][j] =  0.5*fw[1][i][j]; //Dw is zero
			ae[1][i][j] =  max(-fe[1][i][j],D - 0.5*fe[1][i][j]);
			as[1][i][j] = max(fs[1][i][j], D + 0.5*fs[1][i][j]);
			an[1][i][j] = max(-fn[1][i][j], D - 0.5*fn[1][i][j]);		
			sp[1][i][j] = -2.0*D;
			
			a[1][i][j] = aw[1][i][j] + ae[1][i][j] + as[1][i][j] + an[1][i][j] 
			+ fe[1][i][j] - fw[1][i][j] + fn[1][i][j] - fs[1][i][j] - sp[1][i][j] ;
				
			}
			
			else if (j==n) //right boundary
			{
			
			fw[1][i][j] = 0.5*(u_star[i+1][j-1]+u_star[i][j-1]);
			fe[1][i][j] = 0.0;//0.5*(u_star[i+1][j]+u_star[i][j]);
			fs[1][i][j] = 0.5*(v_star[i][j]+v_star[i-1][j]);
			fn[1][i][j] = 0.5*(v_star[i+1][j]+v_star[i][j]);
			
		
			
			aw[1][i][j] =  max(fw[1][i][j],D + 0.5*fw[1][i][j]); 
			ae[1][i][j] = 0.5*fe[1][i][j];/// De is zero
			as[1][i][j] = max(fs[1][i][j],D + 0.5*fs[1][i][j]);
			an[1][i][j] =  max(-fn[1][i][j],D - 0.5*fn[1][i][j]);		
			sp[1][i][j] = -2.0*D;
			
			a[1][i][j] = aw[1][i][j] + ae[1][i][j] + as[1][i][j] + an[1][i][j] 
			+ fe[1][i][j] - fw[1][i][j] + fn[1][i][j] - fs[1][i][j] - sp[1][i][j] ;
				
			}
			
			else if (i==1 && (j>1 && j<n) ) //bottom boundary
			{
			
			fw[1][i][j] = 0.5*(u_star[i+1][j-1]+u_star[i][j-1]);
			fe[1][i][j] = 0.5*(u_star[i+1][j]+u_star[i][j]);
			fs[1][i][j] = 0.5*(v_star[i][j]+v_star[i-1][j]);
			fn[1][i][j] = 0.5*(v_star[i+1][j]+v_star[i][j]);
			
		
			
			aw[1][i][j] = max(fw[1][i][j],D + 0.5*fw[1][i][j]);
			ae[1][i][j] =  max(-fe[1][i][j],D - 0.5*fe[1][i][j]);
			as[1][i][j] = max(fs[1][i][j],D + 0.5*fs[1][i][j]);
			an[1][i][j] =  max(-fn[1][i][j],D - 0.5*fn[1][i][j]);		
			
			
			a[1][i][j] = aw[1][i][j] + ae[1][i][j] + as[1][i][j] + an[1][i][j] 
			+ fe[1][i][j] - fw[1][i][j] + fn[1][i][j] - fs[1][i][j] - sp[1][i][j] ;
			}
			
			else if(i==n-1 && (j>1 && j<n)) // top boundary
			{
			fw[1][i][j] = 0.5*(u_star[i+1][j-1]+u_star[i][j-1]);
			fe[1][i][j] = 0.5*(u_star[i+1][j]+u_star[i][j]);
			fs[1][i][j] = 0.5*(v_star[i][j]+v_star[i-1][j]);
			fn[1][i][j] = 0.5*(v_star[i+1][j]+v_star[i][j]);
			
		
			
			aw[1][i][j] = max(fw[1][i][j],D + 0.5*fw[1][i][j]);
			ae[1][i][j] =  max(-fe[1][i][j],D - 0.5*fe[1][i][j]);
			as[1][i][j] = max(fs[1][i][j],D + 0.5*fs[1][i][j]);
			an[1][i][j] =  max(-fn[1][i][j],D - 0.5*fn[1][i][j]);		
			
			
			a[1][i][j] = aw[1][i][j] + ae[1][i][j] + as[1][i][j] + an[1][i][j] 
			+ fe[1][i][j] - fw[1][i][j] + fn[1][i][j] - fs[1][i][j] - sp[1][i][j] ;
			}
			
			
		}
		
		
	}
	
	// y- momentum equations interior points
		for(int i=2;i<=n-2;i++)
	{
		for(int j=2;j<=n-1;j++)
		{
		
			
			
			fw[1][i][j] = 0.5*(u_star[i+1][j-1]+u_star[i][j-1]);
			fe[1][i][j] = 0.5*(u_star[i+1][j]+u_star[i][j]);
			fs[1][i][j] = 0.5*(v_star[i][j]+v_star[i-1][j]);
			fn[1][i][j] = 0.5*(v_star[i+1][j]+v_star[i][j]);
			
		
			aw[1][i][j] = max(fw[1][i][j],D + 0.5*fw[1][i][j]);
			ae[1][i][j] =  max(-fe[1][i][j],D - 0.5*fe[1][i][j]);
			as[1][i][j] = max(fs[1][i][j],D + 0.5*fs[1][i][j]);
			an[1][i][j] =  max(-fn[1][i][j],D - 0.5*fn[1][i][j]);		
			
			
			a[1][i][j] = aw[1][i][j] + ae[1][i][j] + as[1][i][j] + an[1][i][j] 
			+ fe[1][i][j] - fw[1][i][j] + fn[1][i][j] - fs[1][i][j] - sp[1][i][j] ;
				
		}
	}
}

// function for solving momentum equations
void momentum_eq(int n,double D, double h, double p_star[][n+1], double u_star[n+1][n+1], double v_star[n+1][n+1], double a[2][n+1][n+1], 
		double aw[2][n+1][n+1], double ae[2][n+1][n+1],double as[2][n+1][n+1], double an[2][n+1][n+1],
		 double fw[2][n+1][n+1], double fe[2][n+1][n+1],
				double fs[2][n+1][n+1], double fn[2][n+1][n+1], 
				double sp[2][n+1][n+1], double su[2][n+1][n+1]) 
{

	updating_coefficient(n,D,u_star,v_star,a,aw,ae,as,an,fw,fe,fs,fn,sp,su); // updating momentum equation coefficient
	// U momentum equation 
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n-1;j++)
		{
	
			u_star[i][j] = ((p_star[i][j]-p_star[i][j+1])*h + 
					aw[0][i][j]*u_star[i][j-1] + ae[0][i][j]*u_star[i][j+1] 
					+ as[0][i][j]*u_star[i-1][j] + an[0][i][j]*u_star[i+1][j]+su[0][i][j])/a[0][i][j];
					
			
		}
	}
	
	updating_coefficient(n,D,u_star,v_star,a,aw,ae,as,an,fw,fe,fs,fn,sp,su); // updating momentum equation coefficient
	
		for(int i=1;i<=n-1;i++)
	{
		for(int j=1;j<=n;j++)
		{
		
			
			v_star[i][j] = ((p_star[i][j]-p_star[i+1][j])*h + 
					aw[1][i][j]*v_star[i-1][j] + ae[1][i][j]*v_star[i+1][j] 
					+ as[1][i][j]*v_star[i][j-1] + an[1][i][j]*v_star[i][j+1]+su[1][i][j])/a[1][i][j];
			
					
		}
	}
	
	
	updating_coefficient(n,D,u_star,v_star,a,aw,ae,as,an,fw,fe,fs,fn,sp,su); // updating momentum equation coefficient
}

//function for solving pressure correction equation
void p_dash_equation(int n, double h, double u_star[n+1][n+1], double v_star[n+1][n+1], double p_dash[][n+1], double a[2][n+1][n+1])
{
	double aa,ae,aw,as,an,bij;
	
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
		
			if(i==1 && j==1) //left bottom corner
			{
				aw = 0.0;
				ae = h/a[0][i][j];
				as = 0.0;
				an = h/a[1][i][j];
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
			} 
			
			else if(i==1 && (j>1 && j<n))///bottom boundary
			{
				aw = h/a[0][i][j-1];
				ae = h/a[0][i][j];
				as = 0.0;
				an = h/a[1][i][j];
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j]- v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
			}
			else if(i==1 && j==n) /// bottom right corner
			{
				aw = h/a[0][i][j-1];
				ae = 0.0;
				as = 0.0;
				an = h/a[1][i][j];
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
				
			}
			else if(i==n && (j>1 && j<n)) // top boundary
			{
				aw = h/a[0][i][j-1];
				ae = h/a[0][i][j];
				as = h/a[1][i-1][j];
				an = 0.0;
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
			} 
			
			
			else if(i==n && j==1) //left top corner
			{
				aw = 0.0;
				ae = h/a[0][i][j];
				as = h/a[1][i-1][j];
				an = 0.0;
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
			} 
			
			else if(j==1 && (i>1 && i<n))///left boundary
			{
				aw = 0.0;
				ae = h/a[0][i][j];
				as = h/a[1][i-1][j];
				an = 2.0*h/a[1][i][j];
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j]- v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
			}
			else if(i==n && j==n) /// top right corner
			{
				aw = h/a[0][i][j-1];
				ae = 0.0;
				as = h/a[1][i-1][j];
				an = 0.0;
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
				
			}
			else if( j==n && (i>1 && i<n)) /// right boundary
			{
				aw = h/a[0][i][j-1];
				ae = 0.0;
				as = h/a[1][i-1][j];
				an = h/a[1][i][j];
				bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
				aa = aw+ae+as+an; 
				
				p_dash[i][j] = (aw*p_dash[i][j-1] + ae*p_dash[i][j+1] + as*p_dash[i-1][j] + an*p_dash[i+1][j] + bij)/aa ;
				
			}
			
		}
	}
	
	for(int i=2;i<=n-1;i++)// interior points
	{
		for(int j=2;j<=n-1;j++)
		{
		
			aw = h/a[0][i][j-1];
			ae = h/a[0][i][j];
			as = h/a[1][i-1][j];
			an = h/a[1][i][j];
			bij = (u_star[i][j-1] - u_star[i][j] + v_star[i-1][j] - v_star[i][j]) ; 
			aa = aw+ae+as+an; 
			
			
			p_dash[i][j] = (aw*p_dash[i-1][j] + ae*p_dash[i+1][j] + as*p_dash[i][j-1] + an*p_dash[i][j+1] + bij)/aa ;
		}
	}
	
	
}


//function to evaluate 2nd norm of vector
double rms_termination_condition(int n, double u[n+1][n+1], double uold[n+1][n+1])
{
	double d=0.0;
	
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			
				d=d+pow(u[j][i]-uold[j][i],2);
		}
	}
	
	return sqrt(d/(n+1)*(n+1));
	 
}

//function to evaluate infinite norm
double infinite_norm(int n, double u_star[n+1][n+1], double uold[n+1][n+1])
{
	double temp;
	temp = fabs(u_star[1][1] - uold[1][1]);
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n-1;j++)
		{
			
				if(temp  <  fabs(u_star[i][j] - uold[i][j]))
				{	
					temp = fabs(u_star[i][j] - uold[i][j]);
				}
		}
	}
	return temp;
}


//function to print values on screen
void print_onscreen(int n, double g[n+1][n+1])
{
	printf("\n");
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			
				printf("%0.4f\t col=%d \n",g[i][j],j);
		}
		printf("row=//////////////////////////////////////////////=%d\n",i);
	}
	 
	 printf("\n");
}

double max(double i, double j)
{
	if(i>0 && i>j)
	{
		return i;
	}
	else if (j>0 && j>i)
	{
		return j;
	}
	else return 0;
} 
