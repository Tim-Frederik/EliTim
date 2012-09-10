#ifndef __HALO_H
#define __HALO_H

double int_GSL_integrate_qag2(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter);
double int_GSL_integrate_l(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter);
double int_GSL_integrate_crude(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter);

double alpha_c(const double k1x, const double k1y, const double k2x,const double k2y);
double beta_c(const double k1x, const double k1y, const double k2x,const double k2y);
double fs_2(const double k1x, const double k1y, const double k2x, const double k2y);
double gs_2(const double k1x, const double k1y,const double k2x, const double k2y);
double fs_3(const double k1x, const double k1y, const double k2x, const double k2y, const double k3x, const double k3y);
double p_lin (double k, double a);
double b_lin (const double k1x, const double k1y, const double k2x, const double k2y, const double a);
double t_lin(const double k1x, const double k1y, const double k2x, const double k2y, const double k3x, const double k3y, const double a);
double mstar(double a);
double conc(double m, double a);
double Omega_m_z(double a);
double delta_c(double a);
double delta_vir(double a); 
double rho_vir(double a);
double r_vir(double m, double a);
double r_s(double m, double a);
double radius (double m);
double sigma2_integrand(double k, void * params);
double dnudm(double a, double m);
double sigma2(double m);
double dsigma2dm(double m);
double nu(double m, double a);
double nufnu(double m, double a);
double massfunc(double m, double a);
double bias_norm(double a);
double B1_nonorm (double m,double a);
double B1 (double m, double a);
double int_rho_nfw(double r, void* params);
double u_nfwm(double c,double k, double m, double z);
double u_nfwm_c(double c,double k, double m, double z);
double inner_I0j (double logm, void *para);
double I0j (int j, double k1, double k2, double k3, double k4, double a);
double inner_I1j (double logm, void *para);
double I1j (int j, double k1, double k2, double k3, double a);
double p_1h(double k, double a);
double bi_1h(double k1,double k2, double k3, double a);
double tri_1h(double k1, double k2, double k3, double k4,double a);
double tri_2h (double k1x, double k1y,double k2x, double k2y,double k3x, double k3y,double a); 
/* convergence trispectra */
double project_tri_2h (double l1x, double l1y,double l2x, double l2y,double l3x, double l3y); 
double project_tri_1h (double l1x, double l1y,double l2x, double l2y,double l3x, double l3y); 
double project_tri_lin (double l1x, double l1y,double l2x, double l2y,double l3x, double l3y); 
double project_tri_lin_cov (double l1,double l2, int alpha);
double project_tri_1h_cov (double l1,double l2, int alpha);
double project_tri_2h_cov (double l1,double l2, int alpha);
double HSV_cov(double l1, double l2, double fsky);


double tri_lin_cov(double k1, double k2, double a);
double tri_1h_cov(double k1, double k2, double a);
double tri_2h_cov(double k1, double k2, double a);

double inner_project_tri_1h_cov_tomo(double a, void *params);
double inner_project_tri_lin_cov_tomo(double a,void *params);
double inner_project_tri_2h_cov_tomo(double a,void *params);
double project_tri_2h_cov_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double project_tri_1h_cov_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double project_tri_lin_cov_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double HSV_cov_tomo(double l1, double l2, double fsky, int z1, int z2, int z3, int z4);
double int_for_HSV_tomo(double a, void *params);



#endif
