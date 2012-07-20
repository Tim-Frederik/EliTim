#ifndef __EB_FUNCTIONS_H
#define __EB_FUNCTIONS_H

double integrand_j0_Tfunc(double vartheta, void *params);
double W_calc_j0(double ell,int ORDER,int max,int min);
double integrand_j2_Tfunc(double vartheta, void *params);
double W_calc_j2(double ell,int ORDER,int max,int min);
void Y_minus(double x, double eta, double *ym);
void Y_plus(double x, double eta, double *yp);
double W1(double theta_1, double zeta_1, double zeta_2);
double W2(double theta_2, double zeta_3, double zeta_4);
double integrand2_minus(double theta2, void *params);
double integrand1_minus(double theta1, void *params);
double integrand2_plus(double theta2, void *params);
double integrand1_plus(double theta1, void *params);
double Z_plus(double vartheta, double *array);
double Z_minus(double vartheta,double *array);
void create_Ztable_lin(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z);
void create_Ztable_log(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z);
void create_Tm_COSEBIs_table(double vt_max,double vt_min,char *filename2,int table_Tm_COSEBIs,int N_ORDER);

void normalization(int vartheta_max,int vartheta_min,double *store);
void roots(int vartheta_max,int vartheta_min,double store[20][21]);

double T_lin_plus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min);
double T_lin_minus_integrand_COSEBIS(double t, void *params);
double T_lin_minus_COSEBIS(double theta,int ORDER,double vt_max,double vt_min);
double T_log_plus_COSEBIS(double theta,int ORDER,double vt_max,double vt_min);
double T_log_minus_integrand_COSEBIS(double theta, void *params);
double T_log_minus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min);


double T_plus_MAP(double x);
double T_minus_MAP(double x);
#endif

