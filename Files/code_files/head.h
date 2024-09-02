
double vectorDotproduct(int n,double vec1[n],double vec2[n]);
double residual(int n,int N_j,double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double x[n], double b[n]);
double matrix_vector_product(int n,int N_j,double ap[n], double aw[n], double as[n], double an[n], double ae[n], double d[n], double Ad[n]);
double two_norm(int n, double r[n]);
double ILU_Pinv_r(int n,int N_j, double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double Pinvxr[n]);
double SIP_Pinv_r(int n,int N_j, double ap[n], double aw[n], double as[n], double an[n], double ae[n], double r[n], double Pinvxr[n]);
