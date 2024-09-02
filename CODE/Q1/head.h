double v_cycle (int n, int k);// v-cycle 

double gauss_iteration(int n, int neu, double x[n+1], double b[n+1]); // solving Ax=b equation using gauss siedel iterations for neu iterations

double residual(int n,double r[n+1], double x[n+1], double b[n+1]);

void ristriction(int m, double dh[2*m+1], double d2h[m+1]  );// a is at fine grid , b is at coarse grid

void prolongation(int m, double dh[2*m+1], double d2h[m+1]  ); // a is at fine grid , b is at coarse grid

double two_norm(int n, double r[n+1]); // 2_norm
