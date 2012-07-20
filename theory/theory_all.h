#ifndef __THEORY_H
#define __THEORY_H

#include "fftw3.h"

#define table_N_a     100
#define table_N_k     1000
#define table_N_s    200
#define table_N_klin 200
#define table_N_thetaH 2048
#define table_N_S2    1000
#define table_N_BN    100
#define table_N_SV    100
#define table_N_DS    1000
#define table_N_halo_s  100
#define table_N_xi_nongauss  1000

/* ============================================================ *
 * Cosmology.							*
 * ============================================================ */
void omega_a(double aa,double *om_m,double *om_v);
double omv_vareos(double a);
double growfac(double a, double omm, double omv);
int func_for_growfac(double a,const double y[],double f[],void *params);
double Delta_L(double k);
double Delta_L_tab(double rk);
double p_lin(double k_NL,double a);
double Tsqr_EH_wiggle(double kk);
double sigma_8_sqr_EH_wiggle();
double int_for_sigma_8_EH_wiggle(double k, void *args);
double Tsqr_EH(double k);
double int_for_sigma_8(double k, void *args);
double int_GSL_integrate_qag(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter);
double int_for_wint2_knl(double logk);
double int_for_wint2_neff(double logk);
double int_for_wint2_ncur(double logk);
void wint2(double r,double *sig,double *d1,double *d2, double amp, int onlysig);
void nonlinscale(double,double *,double *,double *,int *);
void halofit(double rk, double rn, double rncur, double rknl, double plin,double om_m, double om_v,double *pnl,double aa);
double slope_NL(double rn, double rncur, double om_m, double om_v);
void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da);
double Delta_nl_halo(double a,double k_NL);
double sigma_8_sqr();
void emu(double *xstar, double *kstar, double *pdelta, int *outtype);
void getH0fromCMB(double *xstar, double *stuff);

// double P_L(double, double);
// void Delta_nl_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da);
// double g(double a);
// double dlog(double x);
// 
// 
// double sigma_8_sqr();
//
/* ============================================================ *
 * P_kappa							*
 * ============================================================ */
double int_for_chi(double a);
double chi(double a);


double g_source(double a);
double f_K(double w);
double G(double a);
double int_for_g(double);
double g_source(double a);
double prob(double z);
double int_for_zdistr(double z);
double int_for_zdistr_mock(double z);
double int_for_p_shear_shear(double a);
double P_shear_shear(double s);
double pf(double z);
double int_for_p_shear_position(double a);
double P_shear_position(double s);
double int_for_p_position_position(double a);
double P_position_position(double s);

/* ============================================================ *
 * correlation function							*
 * ============================================================ */
double xi_shear_shear(int pm, double theta);
double xi_shear_magnification(int pm, double theta);
double xi_shear_position(int pm, double theta);
double xi_magnification_position(int pm, double theta);
double xi_position_position(int pm, double theta);
void xi_via_hankel_shear_shear(double **xi, double *logthetamin, double *logthetamax);
void xi_via_hankel_shear_magnification(double **xi, double *logthetamin, double *logthetamax);
void xi_via_hankel_shear_position(double **xi, double *logthetamin, double *logthetamax);
void xi_via_hankel_magnification_position(double **xi, double *logthetamin, double *logthetamax);
void xi_via_hankel_position_position(double **xi, double *logthetamin, double *logthetamax);

/* ============================================================ *
 * IA							*
 * ============================================================ */
double int_for_p_GI(double a);
double P_GI(double s); 
double xi_GI(int pm, double theta);
void xi_via_hankel_GI(double **xi, double *logthetamin, double *logthetamax);
double int_for_p_II(double a);
double P_II(double s)  ;
double xi_II(int pm, double theta);
void xi_via_hankel_II(double **xi, double *logthetamin, double *logthetamax);


/* ============================================================ *
 * tomography							*
 * ============================================================ */
double int_for_g2(double aprime);

void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc);

#endif
