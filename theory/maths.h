/* ============================================================ *
maths. h								*
 * TE 25.9.2007							*
 * ============================================================ */

#ifndef __MATHS_H
#define __MATHS_H

#include "fftw3.h"

typedef struct {
     double pi, pi_sqr, twopi, ln2, arcmin,lightspeed;
}con;

typedef struct {
     double eps1;
     double eps2;
     double eps3;
}pre;

typedef struct {
     double a_min;
     double Delta_L_klin_min;
     double Delta_L_klin_max;
     double Pdelta_halo_k_min;
     double Pdelta_halo_k_max;
     double P_2_s_min;
     double P_2_s_max;
     double xi_via_hankel_theta_min;
     double xi_via_hankel_theta_max;
     double M_min;
     double M_max;
     double halo_s_min;
     double halo_s_max;
}lim;



static double darg __attribute__((unused)),maxarg1 __attribute__((unused)), maxarg2 __attribute__((unused));

/*DSQR nicht als Produkt schreiben !!! NOT DSQR*DSQR....gives undefined warning */
#define DSQR(a) ( (darg=(a))==0.0 ? 0.0 : darg*darg )  
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) < (maxarg2) ? (maxarg1) : (maxarg2))

#define NR_END 1
#define FREE_ARG char*


 void sm2_error(char *s);

                             
/* ============================================================ *
 * Functions from Numerical Recipes.				*
 * ============================================================ */
double *sm2_vector(long nl, long nh);
int *int_vector(long nl, long nh);
long *long_vector(long nl, long nh);
double **sm2_matrix(long nrl, long nrh, long ncl, long nch);
double sm2_dfridr(double (*func)(double,double), double x, double h, double *err,double aa);
void sm2_spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void sm2_splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double sm2_trapzd(double (*func)(double), double a, double b, int n, double *s);
double sm2_interpol(double *f, int n, double a, double b, double dx, double x,double lower, double upper);
double sm2_interpol2d(double **f, int nx, double ax, double bx, double dx, double x,int ny, double ay, double by, double dy, double y,double lower, double upper);
double sm2_interpol_coyote(double **f,int nx, double ax, double bx, double dx, double x,int ny, double ay, double by, double dy, double y);
double sm2_gammln(double xx);
double sm2_qromb(double (*func)(double), double a, double b);
double sm2_qromo(double (*func)(double), double a, double b,double (*choose)(double(*)(double), double, double, int));
double sm2_midpnt(double (*func)(double), double a, double b, int n);
void sm2_free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void sm2_free_vector(double *v, long nl, long nh);
double bess(double);
double bessj(int ,double);
double bessj1(double);
double bessj0(double);
void cdgamma(fftw_complex x, fftw_complex *res);
void delta_theta_bin_mid(int Ndata,double *theta_binmid,double *deltatheta, int Lin);
double sm2_interpol2d_noextra(double **f,int nx, double ax, double bx, double dx, double x,    int ny, double ay, double by, double dy, double y, double lower, double upper);
#endif
