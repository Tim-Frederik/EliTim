#ifndef __THEORY_H
#define __THEORY_H

#include "fftw3.h"

#define table_N_a     200
#define table_N_k     100
#define table_N_s     100
#define table_N_klin 2000
#define table_N_thetaH 2048

/* ============================================================ *
 * Cosmology.							*
 * ============================================================ */

double Tsqr_EH(double k);
double omv_vareos(double a);
double growfac(double a, double omm, double omv);
int func_for_growfac(double a,const double y[],double f[],void *params);
void omega_a(double aa,double *om_m,double *om_v);
double sigma_R_sqr(double);
double P_L(double, double);
void Delta_nl_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da);
double Delta_nl_halo(double, double);
void nonlinscale(double,double *,double *,double *,int *);
double g(double a);
double Delta_L_tab(double rk);
double Delta_L(double k);
double dlog(double x);
void halofit(double rk, double rn, double rncur, double rknl, double plin,double om_m, double om_v,double *pnl,double aa);
void wint2(double r,double *sig,double *d1,double *d2, double amp, int onlysig);
double int_for_sigma_8(double k, void *args);
double int_GSL_integrate_qag(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter);
double sigma_8_sqr();
void emu(double *xstar, double *kstar, double *pdelta, int *outtype);
void getH0fromCMB(double *xstar, double *stuff);

/* ============================================================ *
 * P_kappa							*
 * ============================================================ */

double g_source(double a);
double f_K(double w);
double G(double a);
double int_for_w(double);
double w(double a);
double prob(double z);
double int_for_zdistr(double z);
double zdistr_bg(double z);
double P_shear_shear(double);
double P_shear_position(double);
double P_position_position(double);
double P2_tomo(double s, double max, double min, double alpha,double beta, double z0);

/* ============================================================ *
 * correlation function							*
 * ============================================================ */
double xi_tomo(int pm, double theta, double max, double min, double alpha,double beta, double z0);
double xi(int pm, double theta);
void xi_via_hankel(double **xi, double *logthetamin, double *logthetamax);
void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc);

#endif
