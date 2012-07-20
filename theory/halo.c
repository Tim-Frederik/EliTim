#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include "theory_all.h"
#include "survey.h"
#include "maths.h"
#include "cosmology.h"
#include "halo.h"
#include "fftw3.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>

extern cosmopara cosmology;
extern sur survey;
extern con constants;
extern lim limits;
extern redshiftshearpara redshiftshear;
extern redshiftclusteringpara redshiftclustering;

/*==============================================================*/
double int_GSL_integrate_qag2(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
 gsl_integration_cquad(&F,a,b,0,1e-5,w,&res,&err,0);
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}

double int_GSL_integrate_l(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
 gsl_integration_cquad(&F,a,b,0,1.e-3,w,&res,&err,0);
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}

double int_GSL_integrate_crude(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *wcrude = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
 gsl_integration_cquad(&F,a,b,0,1.e-2,wcrude,&res,&err,0);
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(wcrude);
  return res;
}

/* ************************Begin Routines for halo model ************* */
/* fundamental mode coupling functions */

double alpha_c(const double k1x, const double k1y,const double k2x,const double k2y)    {
  
  double k1 = sqrt(k1x*k1x+k1y*k1y);
  //        assert(fabs(k1)>1e-10); 
  
  if (fabs(k1) <= 1e-10) return 1e10;             
  return 1.+(k1x*k2x+k1y*k2y)/(k1*k1);
}

double beta_c(const double k1x, const double k1y,const double k2x,const double k2y)
{
  double k1 = sqrt(k1x*k1x+k1y*k1y);
  double k2 = sqrt(k2x*k2x+k2y*k2y+1.e-5);
  double k1k2 = (k1x*k2x+k1y*k2y+1.e-5);
  //        assert((fabs(k1)>1e-10) || (fabs(k2)>1e-10));
  
  const double beta = k1k2*(k1*k1+k2*k2+2.*k1k2)/(2.*k1*k1*k2*k2);
  
  if ((fabs(k1) <= 1e-10) || (fabs(k2) <= 1e-10)) return 1e10;
  return beta;
}

/* F & G Kernels */


double fs_2(const double k1x, const double k1y, const double k2x, const double k2y)
{
  double k1,k2,res;
  k1 = sqrt(k1x*k1x+k1y*k1y+1.e-10);
  k2 = sqrt(k2x*k2x+k2y*k2y+1.e-10);
  res = 5./7.+2./7.*pow((k1x*k2x+k1y*k2y)/(k1*k2),2.)+.5*((k1x*k2x+k1y*k2y)/(k1*k2))*(k1/k2+k2/k1);
  return res;
  
}

double gs_2(const double k1x, const double k1y,const double k2x, const double k2y)
{
  double k1,k2,res;
  k1 = sqrt(k1x*k1x+k1y*k1y+1.e-10);
  k2 = sqrt(k2x*k2x+k2y*k2y+1.e-10);
  res = 3./7.+4./7.*pow((k1x*k2x+k1y*k2y)/(k1*k2),2.)+.5*((k1x*k2x+k1y*k2y)/(k1*k2))*(k1/k2+k2/k1);
  return res;
}

double fs_3(const double k1x, const double k1y, const double k2x, const double k2y, const double k3x, const double k3y)
{
//  printf("fs_31\n");
  double term1,term2,term3;
  term1=7./54.*(alpha_c(k1x,k1y,k2x+k3x,k2y+k3y)*fs_2(k2x,k2y,k3x,k3y) 
  + alpha_c(k2x,k2y,k1x+k3x,k1y+k3y)*fs_2(k1x,k1y,k3x,k3y)
  + alpha_c(k3x,k3y,k1x+k2x,k1y+k2y)*fs_2(k1x,k1y,k2x,k2y));
//  printf("fs_32\n");
  
  term2=7./54.*(alpha_c(k1x+k2x,k1y+k2y,k3x,k3y)*gs_2(k1x,k1y,k2x,k2y) 
  +     alpha_c(k1x+k3x,k1y+k3y,k2x,k2y)*gs_2(k1x,k1y,k3x,k3y)
  +     alpha_c(k2x+k3x,k2y+k3y,k1x,k1y)*gs_2(k2x,k2y,k3x,k3y));
//  printf("fs_33\n");
  
  term3=4./54.*(beta_c(k1x,k1y,k2x+k3x,k2y+k3y)*gs_2(k2x,k2y,k3x,k3y) 
  + beta_c(k2x,k2y,k1x+k3x,k1y+k3y)*gs_2(k1x,k1y,k3x,k3y)
  + beta_c(k3x,k3y,k1x+k2x,k1y+k2y)*gs_2(k1x,k1y,k2x,k2y));
  return term1+term2+term3;
}

double b_lin (const double k1x, const double k1y, const double k2x, const double k2y, const double a)
{
  double k1,k2,k3,k3x,k3y;
  k3x = -k1x-k2x;
  k3y = -k1y-k2y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  return 2.*(fs_2(k1x,k1y,k2x,k2y)*p_lin(k1,a)*p_lin(k2,a)+fs_2(k2x,k2y,k3x,k3y)*p_lin(k2,a)*p_lin(k3,a)+fs_2(k3x,k3y,k1x,k1y)*p_lin(k3,a)*p_lin(k1,a));
}

double t_lin(const double k1x, const double k1y, const double k2x, const double k2y, const double k3x, const double k3y, const double a)
{
 
//  printf("tlin %le %le %le %le %le %le %le\n",k1x,k1y,k2x,k2y,k3x,k3y,a);
  double k1,k2,k3,k4,k4x,k4y;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  double k12,k13,k14;
  
  k12 =sqrt(k1*k1+k2*k2+2.*k1x*k2x+2.*k1y*k2y);
  k13 =sqrt(k1*k1+k3*k3+2.*k1x*k3x+2.*k1y*k3y);
  k14 =sqrt(k1*k1+k4*k4+2.*k1x*k4x+2.*k1y*k4y);
  //         k23 =sqrt(k2*k2+k3*k3+2.*k2x*k3x+2.*k2y*k3y);
  //         k24 =sqrt(k2*k2+k4*k4+2.*k2x*k4x+2.*k2y*k4y);
  //         k34 =sqrt(k3*k3+k4*k4+2.*k3x*k4x+2.*k3y*k4y);
   if(isnan(k12)) k12=0.0;
   if(isnan(k13)) k13=0.0;
   if(isnan(k14)) k14=0.0;

//   printf("%le %le %le %le\n",k12,k13,k14,a);
//   printf("%le\n",p_lin(k12,a));
//   printf("%le\n",p_lin(k13,a));
//   printf("%le\n",p_lin(k14,a));
  
  double fs3_term,fs2_term;
  fs3_term=  fs_3(k1x,k1y,k2x,k2y,k3x,k3y)*p_lin(k1,a)*p_lin(k2,a)*p_lin(k3,a)
  +	fs_3(k1x,k1y,k2x,k2y,k4x,k4y)*p_lin(k1,a)*p_lin(k2,a)*p_lin(k4,a)
  +	fs_3(k1x,k1y,k3x,k3y,k4x,k4y)*p_lin(k1,a)*p_lin(k3,a)*p_lin(k4,a)
  +	fs_3(k2x,k2y,k3x,k3y,k4x,k4y)*p_lin(k2,a)*p_lin(k3,a)*p_lin(k4,a);
//  printf("after fs3\n");
  
  fs2_term=  fs_2(-k1x,-k1y,(k1x+k2x),(k1y+k2y))*fs_2(k3x,k3y,(k1x+k2x),(k1y+k2y))*p_lin(k1,a)*p_lin(k12,a)*p_lin(k3,a)
  +  	fs_2(-k1x,-k1y,(k1x+k2x),(k1y+k2y))*fs_2(k4x,k4y,(k1x+k2x),(k1y+k2y))*p_lin(k1,a)*p_lin(k12,a)*p_lin(k4,a)
  +  	fs_2(-k2x,-k2y,(k1x+k2x),(k1y+k2y))*fs_2(k3x,k3y,(k1x+k2x),(k1y+k2y))*p_lin(k2,a)*p_lin(k12,a)*p_lin(k3,a)
  +  	fs_2(-k2x,-k2y,(k1x+k2x),(k1y+k2y))*fs_2(k4x,k4y,(k1x+k2x),(k1y+k2y))*p_lin(k2,a)*p_lin(k12,a)*p_lin(k4,a)
  +  	fs_2(-k1x,-k1y,(k1x+k3x),(k1y+k3y))*fs_2(k2x,k2y,(k1x+k3x),(k1y+k3y))*p_lin(k1,a)*p_lin(k13,a)*p_lin(k2,a)
  +  	fs_2(-k1x,-k1y,(k1x+k3x),(k1y+k3y))*fs_2(k4x,k4y,(k1x+k3x),(k1y+k3y))*p_lin(k1,a)*p_lin(k13,a)*p_lin(k4,a)
  +  	fs_2(-k3x,-k3y,(k1x+k3x),(k1y+k3y))*fs_2(k2x,k2y,(k1x+k3x),(k1y+k3y))*p_lin(k3,a)*p_lin(k13,a)*p_lin(k2,a)
  +  	fs_2(-k3x,-k3y,(k1x+k3x),(k1y+k3y))*fs_2(k4x,k4y,(k1x+k3x),(k1y+k3y))*p_lin(k3,a)*p_lin(k13,a)*p_lin(k4,a)
  +  	fs_2(-k1x,-k1y,(k1x+k4x),(k1y+k4y))*fs_2(k2x,k2y,(k1x+k4x),(k1y+k4y))*p_lin(k1,a)*p_lin(k14,a)*p_lin(k2,a)
  +  	fs_2(-k1x,-k1y,(k1x+k4x),(k1y+k4y))*fs_2(k3x,k3y,(k1x+k4x),(k1y+k4y))*p_lin(k1,a)*p_lin(k14,a)*p_lin(k3,a)
  +  	fs_2(-k4x,-k4y,(k1x+k4x),(k1y+k4y))*fs_2(k2x,k2y,(k1x+k4x),(k1y+k4y))*p_lin(k4,a)*p_lin(k14,a)*p_lin(k2,a)
  +  	fs_2(-k4x,-k4y,(k1x+k4x),(k1y+k4y))*fs_2(k3x,k3y,(k1x+k4x),(k1y+k4y))*p_lin(k4,a)*p_lin(k14,a)*p_lin(k3,a);
// printf("after fs2\n");
  return 6.*fs3_term+4.*fs2_term;
}

/************** begin halo properties *****************/
/*Bullock et al. concentration prescription - feel free to experiment with more modern fit formulae */
/*******note that mstar is only evaluated at a= 1 so far - redshift dependence is taken care of through Bullock et al. concentration prescription ***/
double mstar(double a)
{
  static double OMEGA_M = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0  = -42.;
  static double SIGMA_8 = -42.;
  
  static double Mstar0;
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec)
  {
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    double  m, s,dc;
    dc = delta_c(1.0);
    m = limits.M_min;
    s = sigma2(m);
    while(s > dc){
      m = m*1.01;
      s = sqrt(sigma2(m));
    }
    printf("calculated mstar(z=0) %e\n", m);
    Mstar0 = m;
  }
  return Mstar0;
}
/******** use c_* = 10, alpha = 0.2, as proposed in Takada & Jain 2003 **********/
double conc(double m, double a)
{
  return 10.*pow(mstar(1.0)/m,.2);
}
/* routines to convert scales/radii to halo masses, and vice versa */
double Omega_m_z(double a)
{
  return cosmology.Omega_m*pow(a,-3.0)/ (cosmology.Omega_m*pow(a,-3.0)+cosmology.Omega_v);
}

double delta_c(double a) /*see Henry 2000 eq. A16 == Nakamura & Suto eq. C28 */
{
  double x= a*a*a*(1./cosmology.Omega_m -1);
  return 3.0/20.0*pow(12.0*M_PI,2./3.)*(1.0-0.0123*log(1+x));
}
double delta_vir(double a)                 /* see Nakamura & Suto 1997 eq. C19 == Henry 2000 eq. A17 */
{
  double x= a*a*a*(1./cosmology.Omega_m -1);
  return 18.*M_PI*M_PI * (1.+0.4093*pow(x,0.90524));  
}

double rho_vir(double a) //virial density in solar masses/h/(H0/c)^3
{
  
  return (delta_vir(a) *cosmology.rho_m);
}

double r_vir(double m, double a) //calculate r_vir in c/H0 given m in (solar masses/h)
{
  return pow(3./(4.*M_PI)*(m/rho_vir(a)),.333333333); 
}

double r_s(double m, double a)
{
  return r_vir(m,a)/conc(m,a);
}

double radius(double m)
{
  double r = pow(3./4.*m/(M_PI*cosmology.rho_m),.333333333);
  return r;
}

/*++++++++++++++++++++++++++++++++++++++++*
 *  Variance of the density field         *
 *++++++++++++++++++++++++++++++++++++++++*/


double sigma2_integrand(double k, void * params)   // inner integral
{
  double *array = (double*)params;
  double m = array[0];
  double x= k*radius(m);
  
  double inner= p_lin(k,1.0)*k*k*pow(3./(x*x*x) * (sin(x) - x*cos(x)),2.)/(2.*constants.pi*constants.pi);
  //printf("%le %le\n",k,p_lin(k,1.0));
  return inner;
}

double dnudm(double a, double m){
  return log(nu(exp(m),a));
}

double sigma2(double m)
{
  static double OMEGA_M = -42.;
  static double N_SPEC  = -42.;
  static double SIGMA_8 = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  
  static double table_S2[table_N_S2];
  static double dm = .0, logmmin =1., logmmax = 1.;
  double mlog,f1,result,array[1];
  int i;
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec)
  {
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    logmmin = log(limits.M_min/10.0);
    logmmax = log(limits.M_max*10.0);
    dm = (logmmax - logmmin)/(table_N_S2-1.);
    mlog = logmmin;
    
    for (i=0; i<table_N_S2; i++, mlog += dm) {  
      array[0] = exp(mlog);
      result = int_GSL_integrate_qag2(sigma2_integrand,(void*)array,1.e-6,1.0e+6,NULL,5000);
      table_S2[i]=log(result);
    }
  }
  
  mlog=log(m);
  f1=sm2_interpol(table_S2, table_N_S2, logmmin, logmmax, dm,mlog, 1.0,1.0 );
  return exp(f1);  
}

double dsigma2dm(double m)
{     
  static double OMEGA_M = -42.;
  static double N_SPEC  = -42.;
  static double SIGMA_8 = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  
  static double table_DS[table_N_DS];
  static double dm = .0, logmmin = 1.0,logmmax = 1.0;
  double mlog,f1,result,error;
  int i;
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec)
  {
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    logmmin = log(limits.M_min/10.0);
    logmmax = log(limits.M_max*10.0);
    
    dm = (logmmax - logmmin)/(table_N_S2-1.);
    mlog = logmmin;
    
    for (i=0; i<table_N_S2; i++, mlog += dm) {  
      result = sm2_dfridr(dnudm, mlog, 0.01*mlog,&error, 0.9999);
      table_DS[i]=log(result);
    }
  }
  mlog=log(m);
  f1=sm2_interpol(table_DS, table_N_DS, logmmin, logmmax, dm,mlog, 1.0,1.0 );
  return exp(f1);  
}
/*++++++++++++++++++++++++++++++++++++++++*
 *  Mass function (Sheth & Tormen)        *
 *++++++++++++++++++++++++++++++++++++++++*/

double nu(double m, double a)   {
  const double growth_z = growfac(a,1.,1.)/growfac(1.,1.,1.);
  double res = delta_c(a)*delta_c(a)/(sigma2(m)*pow(growth_z,2.));
  return res;
}

double nufnu(double m, double a)
{
  const double q_nu_m_z = 0.707*nu(m,a);
  //return 0.322*(1.+pow(q_nu_m_z,-0.3)) * sqrt(2.0*q_nu_m_z/M_PI) * exp(-q_nu_m_z/2.);
  return 0.19*(1.+pow(q_nu_m_z,-0.3)) * sqrt(2.0*q_nu_m_z/M_PI) * exp(-q_nu_m_z/2.);
}

double massfunc(double m, double a)
{
  double massfct = nufnu(m,a)*cosmology.rho_m/m/m*dsigma2dm(m)/2.0; //dsigma2dm(m) = 2 dlog(nu)/dlogm
  return massfct;
}

//b_1 - normalization due to mass cuts (to get large scale clustering close to PT result)
double bias_norm_integrand(double lnm, void * params)   // inner integral
{
  double *array = (double*)params;
  double a = array[0], m = exp(lnm);
  return massfunc(m,a)/cosmology.rho_m * m*m * B1_nonorm(m,a);
}
double bias_norm(double a)
{
  static double OMEGA_M = -42.;
  static double N_SPEC  = -42.;
  static double SIGMA_8 = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  
  static double table_BN[table_N_BN];
  static double da = .01, amin = 0.1, amax =1.0;
  double aa,result,array[1];
  int i;
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec)
  {
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    
    da = (amax - amin)/(table_N_BN-1.);
    aa= amin;
    
    for (i=0; i<table_N_BN; i++, aa += da) {  
      array[0] = aa;
      result = int_GSL_integrate_qag2(bias_norm_integrand,(void*)array,log(limits.M_min),log(limits.M_max),NULL,5000);
      table_BN[i]=result;
    }
  }
  return sm2_interpol(table_BN, table_N_BN, amin, amax, da,a, 1.0,1.0 );
  
}


//b_1, Sheth& Tormen:  the following routines need to be identical, except for bias_norm factor
double B1_nonorm (double m,double a){
  double q = 0.707,p =0.3; 
  double res = 1. +  (q*nu(m,a) - 1.)/delta_c(a)+2.*p/delta_c(a)/(1.+pow(q*nu(m,a),p));
  return res;
}

double B1 (double m,double a){
  double q = 0.707,p =0.3; 
  double res = 1. +  (q*nu(m,a) - 1.)/delta_c(a)+2.*p/delta_c(a)/(1.+pow(q*nu(m,a),p));
  return res;/*/bias_norm(a);*/
}
/*++++++++++++++++++++++++++++++++++++++++*
 *  FT of NFW Profile (CS 81)             *
 *++++++++++++++++++++++++++++++++++++++++*/

double int_rho_nfw(double r, void *params){
  double *array = (double*)params;
  double k = array[0];
  double c = array[3];
  double rv = array[4];
  double rs =rv/c;
  /* return Fourier kernel * NFW profile */
  return r*sinl(k*r)/k*pow(rs,-3.)*1./(log(1.+c)-c/(1.+c))/(r/rs)*pow(1+r/rs,-2.)*1.e+5;
}

double u_nfw(double c, double k, double m, double a){
  // FFT density of NFW profile, truncated at r_vir, through direct integration
  // use this routine in inner I[0-1]j and modify int_rho_nfw if you want to try different halo profiles 
  // e.g. to include adiabatic contraction, AGN feedback, etc.
  
  double array[5],val;
  array[0] = k;
  array[1] = m;
  array[2] = a;
  array[3] = c;
  array[4] = r_vir(m,a);
  
  val = int_GSL_integrate_qag2(int_rho_nfw, (void*)array, 0,array[4],NULL, 5000);
  return val;
}


double u_nfw_c(double c,double k, double m, double aa)
{
  // analytic FT of NFW profile, from Cooray & Sheth 01
  const double x = k * r_vir(m,aa)/c;
  double xl, xu;
  xl = x; xu = (1.+c)*x;
  const double val = (sinl(xl)*(gsl_sf_Si(xu)-gsl_sf_Si(xl))- sinl(c*x)/((1.+c)*x) +cosl(xl)*(gsl_sf_Ci(xu)-gsl_sf_Ci(xl)))*1./(log(1.+c)-c/(1.+c)) *1.e+5;
  return val;
}


double inner_I0j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
//printf("%le\n",conc(m,a));
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
//  printf("%le %le %le %le\n",array[0],array[1],array[2],array[3]);
 for (l = 0; l< j; l++){
    u = u*u_nfw_c(c,array[l],m,a);
  }
 
  return massfunc(m,a)*m*pow(m/cosmology.rho_m,(double)j)*u;
}

double I0j (int j, double k1, double k2, double k3, double k4,double a){
  double array[7];
  double result,error;
  gsl_integration_cquad_workspace *w;
  gsl_function H;
    
  array[0]=k1;
  array[1]=k2;
  array[2]=k3;
  array[3]=k4;
  array[4]=0;
  array[5]=(double)j;
  array[6]=a;
  
  w = gsl_integration_cquad_workspace_alloc (2000);
  
  H.function = &inner_I0j;
  H.params = (void*)array;
  
 gsl_integration_cquad (&H,log(limits.M_min),log(limits.M_max), 0, 1.e-3, w,&result,&error,0);  // changed from 1e1 and 1e20
  gsl_integration_cquad_workspace_free(w);
  return result;
}

double inner_I1j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
  for (l = 0; l< j; l++){
    u = u*u_nfw_c(c,array[l],m,a);
  }
  return massfunc(m,a)*m*pow(m/cosmology.rho_m,(double)j)*u*B1(m,a);
}

double I1j (int j, double k1, double k2, double k3,double a){
  double array[7];
  double result,error;
  gsl_integration_cquad_workspace *w;
  gsl_function H;
  
  
  array[0]=k1;array[1]=k2;array[2]=k3;
  array[4]=0;array[5]=(double)j,
  array[6]=a;
  
  w = gsl_integration_cquad_workspace_alloc (2000);
  
  H.function = &inner_I1j;
  H.params = (void*)array;
  
 gsl_integration_cquad (&H,log(limits.M_min),log(limits.M_max), 0, 1.e-3, w, &result, &error,0);  // changed from 1e1 and 1e20
  gsl_integration_cquad_workspace_free(w);
  return result;
}

/*++++++++++++++++++++++++++++++++++++++++*
 * 1Halo Terms                            *
 *++++++++++++++++++++++++++++++++++++++++*/

double p_1h(double k, double a)
{
  return I0j(2,k,k,0,0,a)*1.e-10; 
}
double bi_1h(double k1,double k2, double k3, double a)
{
  return I0j(3,k1,k2,k3,0,a)*1.e-15;
}
double tri_1h(double k1, double k2, double k3, double k4,double a)
{
  return I0j(4,k1,k2,k3,k4,a)*1.e-20;         
}

/*++++++++++++++++++++++++++++++++++++++++*
 * 2Halo Terms                            					   *
 *++++++++++++++++++++++++++++++++++++++++*/
double p_2h(double k, double a)
{
  return pow(I1j(1,k,0,0,a),2.0)*p_lin(k,a)*1.e-10; 
}
double tri_2h_13 (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  double k1,k2,k3,k4,k4x,k4y;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  return (I1j(3,k1,k2,k3,a)*I1j(1,k4,0,0,a)*p_lin(k4,a) + I1j(3,k1,k2,k4,a)*I1j(1,k3,0,0,a)*p_lin(k3,a) + I1j(3,k1,k3,k4,a)*I1j(1,k2,0,0,a)*p_lin(k2,a) + I1j(3,k4,k2,k3,a)*I1j(1,k1,0,0,a)*p_lin(k1,a))*1.e-20;
}

double tri_2h_22 (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  double k1,k2,k3,k4,k4x,k4y,k12,k13,k14;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  k12 = sqrt(k1*k1 + 2.*k1x*k2x+2.*k1y*k2y+k2*k2);
  k13 = sqrt(k1*k1 + 2.*k1x*k3x+2.*k1y*k3y+k3*k3);
  k14 = sqrt(k1*k1 + 2.*k1x*k4x+2.*k1y*k4y+k4*k4);
  
  if(isnan(k12)) k12=0.0;
  if(isnan(k13)) k13=0.0;
  if(isnan(k14)) k14=0.0;

  return (I1j(2,k1,k2,0,a)*I1j(2,k3,k4,0,a)*p_lin(k12,a) + I1j(2,k1,k3,0,a)*I1j(2,k2,k4,0,a)*p_lin(k13,a) + I1j(2,k1,k4,0,a)*I1j(2,k3,k2,0,a)*p_lin(k14,a))*1.e-20;
}

double tri_2h (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  return tri_2h_22 (k1x, k1y, k2x, k2y, k3x, k3y, a) + tri_2h_13 (k1x, k1y, k2x, k2y, k3x, k3y, a);
}


/* ************************* Halo model convergence trispectra ******* */
// static double inner_project_tomo_tri_lin (double a,void *params)
// {
//   
//   double k1x,k1y,k2x,k2y,k3x,k3y,hoverh0,res,weight,wa;
//   double *ar = (double *) params;
//   wa = chi(a);
//   k1x = ar[0]/wa;
//   k1y = ar[1]/wa;
//   k2x = ar[2]/wa;
//   k2y = ar[3]/wa;
//   k3x = ar[4]/wa;
//   k3y = ar[5]/wa;
//   const double m04 = t_lin(k1x,k1y,k2x,k2y,k3x,k3y,a);
//   weight = 9./4.*cosmology.Omega_m*cosmology.Omega_m/(a*a)*(1-wa/w(ar[6]))*3./2.*cosmology.Omega_m/a*(1-wa/w(ar[7]));
//   weight = pow(weight,2.);
//   hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
//   res= weight/pow(wa,2.)/hoverh0;
//   res = res*m04;
//   return res/(a*a); //a^-2 from change of integration variable dchi -> da
// }

static double inner_project_tri_lin (double a,void *params)
{
  
  double k1x,k1y,k2x,k2y,k3x,k3y,hoverh0,res,weight,wa;
  double *ar = (double *) params;
  wa = chi(a);
  k1x = ar[0]/wa;
  k1y = ar[1]/wa;
  k2x = ar[2]/wa;
  k2y = ar[3]/wa;
  k3x = ar[4]/wa;
  k3y = ar[5]/wa;
  //printf("t_lin %le %le %le %le %le %le\n",k1x,k1y,k2x,k2y,k3x,k3y);
  const double m04 = t_lin(k1x,k1y,k2x,k2y,k3x,k3y,a);
  weight = 3./2.*cosmology.Omega_m/a*g_source(a);
  weight = pow(weight,4.);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= weight/pow(wa,2.)/hoverh0;
  res = res*m04;
  return res/(a*a);
}

// static double inner_project_tomo_tri_1h (double a, void *params)
// {
//   double k1,k2,k3,k4,k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,hoverh0,res,weight,wa;
//   double *ar = (double *) params;
//   wa = chi(a);
//   k1x = ar[0]/wa;
//   k1y = ar[1]/wa;
//   k2x = ar[2]/wa;
//   k2y = ar[3]/wa;
//   k3x = ar[4]/wa;
//   k3y = ar[5]/wa;
//   k4x = -(ar[0]+ar[2]+ar[4])/wa; k4y = -(ar[1]+ar[3]+ar[5])/wa; 
//   k1 = sqrt(k1x*k1x + k1y*k1y);
//   k2 = sqrt(k2x*k2x+k2y*k2y);
//   k3 = sqrt(k3x*k3x+k3y*k3y);
//   k4 = sqrt(k4x*k4x+k4y*k4y);
//   const double m04 = tri_1h(k1,k2,k3,k4,a);
//   weight = 9./4.*cosmology.Omega_m*cosmology.Omega_m/(a*a)*(1-wa/w(ar[6]))*3./2.*cosmology.Omega_m/a*(1-wa/w(ar[7]));
//   weight = pow(weight,2.);
//   hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
//   res= weight/pow(wa,2.)/hoverh0;
//   res = res*m04;
//   return res/(a*a);
// }
static double inner_project_tri_1h (double a, void *params)
{
  double k1,k2,k3,k4,k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,hoverh0,res,weight,wa;
  double *ar = (double *) params;
  wa = chi(a);
  k1x = ar[0]/wa;
  k1y = ar[1]/wa;
  k2x = ar[2]/wa;
  k2y = ar[3]/wa;
  k3x = ar[4]/wa;
  k3y = ar[5]/wa;
  k4x = -(ar[0]+ar[2]+ar[4])/wa; k4y = -(ar[1]+ar[3]+ar[5])/wa; 
  k1 = sqrt(k1x*k1x + k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  const double m04 = tri_1h(k1,k2,k3,k4,a);
  weight = 3./2.*cosmology.Omega_m/a*g_source(a);
  weight = pow(weight,4.);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= weight/pow(wa,2.)/hoverh0;
  res = res*m04;
  return res/(a*a);
}

// static double inner_project_tomo_tri_2h (double a,void *params)
// {
//   
//   double k1x,k1y,k2x,k2y,k3x,k3y,hoverh0,res,weight,wa;
//   double *ar = (double *) params;
//   wa = chi(a);
//   k1x = ar[0]/wa;
//   k1y = ar[1]/wa;
//   k2x = ar[2]/wa;
//   k2y = ar[3]/wa;
//   k3x = ar[4]/wa;
//   k3y = ar[5]/wa;
//   const double m04 = tri_2h(k1x,k1y,k2x,k2y,k3x,k3y,a);
//   weight = 9./4.*cosmology.Omega_m*cosmology.Omega_m/(a*a)*(1-wa/w(ar[6]))*3./2.*cosmology.Omega_m/a*(1-wa/w(ar[7]));
//   weight = pow(weight,2.);
//   hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
//   res= weight/pow(wa,2.)/hoverh0;
//   res = res*m04;
//   return res/(a*a);
// }

static double inner_project_tri_2h (double a,void *params)
{
  
  double k1x,k1y,k2x,k2y,k3x,k3y,hoverh0,res,weight,wa;
  double *ar = (double *) params;
  wa = chi(a);
  k1x = ar[0]/wa;
  k1y = ar[1]/wa;
  k2x = ar[2]/wa;
  k2y = ar[3]/wa;
  k3x = ar[4]/wa;
  k3y = ar[5]/wa;
  const double m04 = tri_2h(k1x,k1y,k2x,k2y,k3x,k3y,a);
  weight = 3./2.*cosmology.Omega_m/a*g_source(a);
  weight = pow(weight,4.);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= weight/pow(wa,2.)/hoverh0;
  res = res*m04;
  return res/(a*a);
}

double project_tri_2h (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y){
  double as, array[6],res;
  as = 1./(redshiftshear.z0 +1.);
  array[0] = l1x;
  array[1] = l1y;
  array[2] = l2x;
  array[3] = l2y;
  array[4] = l3x;
  array[5] = l3y;
  res =int_GSL_integrate_crude(inner_project_tri_2h,(void*)array,as,0.9999,NULL,1000);
  return res;
}

double project_tri_1h (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y){
  double as, array[6],res;
  as = 1./(redshiftshear.z0 +1.);
  array[0] = l1x;
  array[1] = l1y;
  array[2] = l2x;
  array[3] = l2y;
  array[4] = l3x;
  array[5] = l3y;
  res =int_GSL_integrate_crude(inner_project_tri_1h,(void*)array,as,0.9999,NULL,1000);
  return res;
}
double project_tri_lin (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y){
  double as, array[6],res;
  as = 1./(redshiftshear.z0 +1.);
  array[0] = l1x;
  array[1] = l1y;
  array[2] = l2x;
  array[3] = l2y;
  array[4] = l3x;
  array[5] = l3y;
  res =int_GSL_integrate_crude(inner_project_tri_lin,(void*)array,as,0.9999,NULL,1000);
  return res;
}
// double project_tomo_tri_2h (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y,int nzi, int nzj){
//   double as, array[8],res;
//   if (nzi <= nzj){as = 1./(zbins.zi[nzj] +1.);}
//   if (nzi > nzj){as = 1./(zbins.zi[nzi] +1.);}
//   array[0] = l1x;
//   array[1] = l1y;
//   array[2] = l2x;
//   array[3] = l2y;
//   array[4] = l3x;
//   array[5] = l3y;
//   array[6] = 1./(zbins.zi[nzi]+1.);
//   array[7] = 1./(zbins.zi[nzj]+1.);
//   res =int_GSL_integrate_crude(inner_project_tomo_tri_2h,(void*)array,as,0.9999,NULL,1000);
//   return res;
// }
// double project_tomo_tri_1h (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y,int nzi, int nzj){
//   double as, array[8],res;
//   if (nzi <= nzj){as = 1./(zbins.zi[nzj] +1.);}
//   if (nzi > nzj){as = 1./(zbins.zi[nzi] +1.);}
//   array[0] = l1x;
//   array[1] = l1y;
//   array[2] = l2x;
//   array[3] = l2y;
//   array[4] = l3x;
//   array[5] = l3y;
//   array[6] = 1./(zbins.zi[nzi]+1.);
//   array[7] = 1./(zbins.zi[nzj]+1.);
//   res =int_GSL_integrate_crude(inner_project_tomo_tri_1h,(void*)array,as,0.9999,NULL,1000);
//   return res;
// }
// double project_tomo_tri_lin (double l1x,double l1y,double l2x,double l2y,double l3x,double l3y,int nzi, int nzj){
//   double as, array[8],res;
//   if (nzi <= nzj){as = 1./(zbins.zi[nzj] +1.);}
//   if (nzi > nzj){as = 1./(zbins.zi[nzi] +1.);}
//   array[0] = l1x;
//   array[1] = l1y;
//   array[2] = l2x;
//   array[3] = l2y;
//   array[4] = l3x;
//   array[5] = l3y;
//   array[6] = 1./(zbins.zi[nzi]+1.);
//   array[7] = 1./(zbins.zi[nzj]+1.);
//   res =int_GSL_integrate_crude(inner_project_tomo_tri_lin,(void*)array,as,0.9999,NULL,1000);
//   return res;
// }

/****************** Halo sample variance term (Sato et al. 2009) **********/

/***convergence of projected NFW profiles (Takada & Jain 2003, Eqs. (26-28)***/
double Gkappa(double x, double c){
  double res =0.0;
  if (x<= c ){
    if ((x-1.0)> 1.e-14){ res = -pow(c*c-x*x,0.5)/((1.0-x*x)*(1.0+c))-1.0/pow(x*x-1.0,1.5)*acos((x*x+c)/(x*(1.0+c)));}
    else if ((1.0-x)> 1.e-14 && x >= 1.e-14){ res = -pow(c*c-x*x,0.5)/((1.0-x*x)*(1.0+c))+1.0/pow(1.0-x*x,1.5)*acosh((x*x+c)/(x*(1.0+c)));}
    else if (x<1.e-14){res = Gkappa(1.e-14,c);}
    else {res =0.5*(Gkappa(0.999,c)+Gkappa(1.001,c));}
  }
  return res;
}
double kappa_nfw(double theta, double m, double a){
  double c, wa,f,rv;
  wa= chi(a);
  c = conc(m,a); 
  rv = r_vir(m,a);
  f = 1./(log(1.0+c)-c/(1.0+c));
  //	return 3./2.*cosmology.Omega_m*g_source(a)*wa/a*m/cosmology.rho_m*f*c*c/(constants.twopi*rv*rv)*Gkappa(c*theta*wa/rv,c);	
  return cosmology.Omega_m*g_source(a)*wa/a*delta_vir(a)*rv*f*c*c*Gkappa(c*theta*wa/rv,c);	
}
double int_for_kappa_l_nfw(double theta, void *params){
  double *ar = (double *) params;
  //return constants.twopi*theta*gsl_sf_bessel_J0(theta*ar[5])* 3./2.*cosmology.Omega_m*g_source(a)*wa/a*m/cosmology.rho_m*f*c*c/(constants.twopi*rv*rv)*Gkappa(c*theta*wa/rv,c);	
  // return constants.twopi*theta*gsl_sf_bessel_J0(theta*ar[5])* cosmology.Omega_m*g_source(a)*wa/a*delta_vir(a)*rv*f*c*c*Gkappa(c*theta*wa/rv,c);	
  return constants.twopi*theta*gsl_sf_bessel_J0(theta*ar[5])*cosmology.Omega_m*g_source(ar[0])*ar[1]/ar[0]*delta_vir(ar[0])*ar[2]*ar[4]*ar[3]*ar[3]*Gkappa(ar[3]*theta*ar[1]/ar[2],ar[3]);	
}
double kappa_l_nfw(double l, double m, double a){
  //	printf("called kappa_l\n");
  double array[6],res =0;
  double d,r,r2,thetav;
  unsigned int n = 1;
  array[0] = a;
  array[1]= chi(a);
  array[2] = r_vir(m,a);
  array[3] =conc(m,a); 
  array[4] = 1./(log(1.0+array[3])-array[3]/(1.0+array[3]));
  array[5] = l;
  thetav =array[2]/array[1];
  r = gsl_sf_bessel_zero_J0(n)/l;
  if (r > thetav){
    res = int_GSL_integrate_qag2(int_for_kappa_l_nfw,(void*)array,1.e-3*r,r,NULL,1000);
    n++;
    r2 = gsl_sf_bessel_zero_J0(n)/l;
    d = res;
    while (fabs(d)> 1.e-3*fabs(res) && r2< thetav){
      d = int_GSL_integrate_qag2(int_for_kappa_l_nfw,(void*)array,r,r2,NULL,1000);
      res +=d;
      r = r2; 
      n++;
      r2 = gsl_sf_bessel_zero_J0(n)/l;
    }
    d = int_GSL_integrate_qag2(int_for_kappa_l_nfw,(void*)array,r,thetav,NULL,1000);
    res +=d;
  } else {res = int_GSL_integrate_qag2(int_for_kappa_l_nfw,(void*)array,1.e-5*thetav,thetav,NULL,1000);}
  //res = int_GSL_integrate_qag2(int_for_kappa_l_nfw,(void*)array,array[2]/array[1]*1.e-4,array[2]/array[1],NULL,1000); //array[2]/array[1] = theta_vir
  //	printf("%e %e %e %e %i\n", m,thetav,l,res,n);
  return res;
  
}

double int_for_kappa_mean (double logm, void *params){
  double *ar = (double *) params;
  double l = ar[0];
  double a = ar[1];
  return exp(logm)*massfunc(exp(logm),a)*B1(exp(logm),a)*pow(kappa_l_nfw(l, exp(logm),a),2.0);
}
double kappa_mean (double l, double a){
  double array[2],res;
  array[0] = l;
  array[1] = a;
  res =int_GSL_integrate_l(int_for_kappa_mean,(void*)array,log(limits.M_min*100.0),log(limits.M_max),NULL,1000);
  return res;
}

double int_for_variance (double logk, void *params){
  double *ar = (double *) params;
  double k = exp(logk);
  double x = pow(4.0*ar[1],0.5)*k*chi(ar[0]); //theta_s*k*chi(a)
  return k*k/constants.twopi*p_lin(k,ar[0])*pow(gsl_sf_bessel_J1(x)/x,2.0);
}
double survey_variance (double a, double fsky){
  static double OMEGA_M = -42.;
  static double OMEGA_B = -42.;
  static double N_SPEC  = -42.;
  static double SIGMA_8 = -42.;
  static double H0 = -42.;
  static double FSKY = -42.;
  
  static double table_SV[table_N_SV];
  static double da = .01, amin = 0.3, amax =0.9999;
  double aa,result,array[2];
  int i;
  if (OMEGA_M != cosmology.Omega_m || OMEGA_B != cosmology.omb ||H0 != cosmology.h0 || FSKY!= fsky|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec )
  {
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    H0=cosmology.h0;
    FSKY=fsky;
    OMEGA_M=cosmology.Omega_m;
    OMEGA_B=cosmology.omb;
    
    da = (amax - amin)/(table_N_SV-1.);
    aa= amin;
    array[1] = fsky;
    
    for (i=0; i<table_N_SV; i++, aa += da) {  
      array[0] = aa;
      result = int_GSL_integrate_qag2(int_for_variance,(void*)array,log(1.e-6),log(1.e+3),NULL,2000);
      table_SV[i]=result;
    }
    
  }
  return sm2_interpol(table_SV, table_N_SV, amin, amax, da,a, 1.0,1.0 );
}

double int_for_HSV(double a, void *params){
  double *ar = (double *) params;
  double res,hoverh0;
  hoverh0 =  sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res = pow(chi(a),4.0)*hoverh0/(a*a); //factor hoverh0/(a*a) from chain rule
  res = res*kappa_mean(ar[0],a)*kappa_mean(ar[1],a);
  res = res*survey_variance(a,ar[2]);
  return res;
}

double HSV_cov(double l1, double l2, double fsky){
  double array[3],res = 0.0;
  array[0] = l1;
  array[1] = l2;
  array[2] = fsky;
  if (sqrt(l1*l2) >= 1.e+2){
    res =int_GSL_integrate_crude(int_for_HSV,(void*)array,1./(1.0+redshiftshear.zdistrpar_zmax),0.99,NULL,1000);
  }
  return res;
}

/******* trispectrum terms for covariance (parallelogram) configurations *******/
double inner_P_bin(double theta, void *params){
  double *ar = (double *) params;
  return p_lin(sqrt(ar[0]*ar[0]+2.0*cos(theta)*ar[0]*ar[1]+ar[1]*ar[1]),ar[2]);
}
double inner_tri_lin_cov(double theta, void *params){
  double res;
  double *ar = (double *) params;
  //printf("%.10le %le %le %le\n",theta,ar[0],ar[1],ar[2]);
  //printf("%le %le\n",cos(theta),sin(theta));
  res=t_lin(ar[0],0, ar[1]*cos(theta), ar[1]*sin(theta), -ar[0],0,ar[2]);
  //printf("%.10le %le\n",theta,res);
  return res;
}
double tri_2h_13_cov (double k1,double k2,double a){
  return (I1j(3,k1,k2,k2,a)*I1j(1,k1,0,0,a)*p_lin(k1,a)+I1j(3,k1,k1,k2,a)*I1j(1,k2,0,0,a)*p_lin(k2,a))*1.e-20;
}

double tri_2h_22_cov (double k1, double k2, double a){
  double array[3];
  array[0] = k1; array[1] = k2; array[2] = a;
  return 2.0/constants.twopi*pow(I1j(2,k1,k2,0,a),2.0)*1.e-20*int_GSL_integrate_crude(inner_P_bin,(void*)array,0,constants.twopi,NULL,1000);
}
double tri_1h_cov(double k1, double k2, double a){
  //printf("%le %le %le\n",k1,k2,a);
  return I0j(4,k1,k1,k2,k2,a)*1.e-20;
}
double tri_2h_cov (double k1,double k2,double a){
  return (tri_2h_22_cov (k1,k2, a)+ tri_2h_13_cov (k1,k2, a));
}
double tri_lin_cov(double k1, double k2, double a){
  double array[3];
  double test;
  array[0] = k1; array[1] = k2;array[2] = a;
  //printf("%le %le %le\n",array[0],array[1],array[2]);
  test= int_GSL_integrate_crude(inner_tri_lin_cov,(void*)array,0,constants.twopi,NULL,1000)/constants.twopi;
  //printf("%le\n",test);
  return test;
}

/* ************************* Halo model convergence trispectra - covariance configuration SHEAR SHEAR******* */
static double inner_project_tri_lin_cov(double a,void *params)
{
  double k1,k2,hoverh0,res,wa,alpha,weight1,weight2,weight3;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  alpha=ar[2];
  const double m04 = tri_lin_cov(k1,k2,a);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  weight1 = pow((g_source(a)*wa/a),alpha);
  weight2 = pow((a*a*hoverh0),(3.0-alpha));
  weight3 = pow((pf(1./a-1.)*cosmology.bias),(4.0-alpha));
  res= weight1*weight2*weight3/pow(wa,6.);
  res = res*m04;
  return res; 
}

static double inner_project_tri_1h_cov(double a, void *params)
{
  double k1,k2,hoverh0,res,wa,alpha,weight1,weight2,weight3;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  alpha=ar[2];
  const double m04 = tri_1h_cov(k1,k2,a);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  weight1 = pow((g_source(a)*wa/a),alpha);
  weight2 = pow((a*a*hoverh0),(3.0-alpha));
  weight3 = pow((pf(1./a-1.)*cosmology.bias),(4.0-alpha));
  res= weight1*weight2*weight3/pow(wa,6.);
  res = res*m04;
  return res; 
}

static double inner_project_tri_2h_cov(double a,void *params)
{
  double k1,k2,hoverh0,res,wa,alpha,weight1,weight2,weight3;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  alpha=ar[2];
  const double m04 = tri_2h_cov(k1,k2,a);
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  weight1 = pow((g_source(a)*wa/a),alpha);
  weight2 = pow((a*a*hoverh0),(3.0-alpha));
  weight3 = pow((pf(1./a-1.)*cosmology.bias),(4.0-alpha));
  res= weight1*weight2*weight3/pow(wa,6.);
  res = res*m04;
  return res; 
}

double project_tri_2h_cov(double l1,double l2,int alpha){
  static double INDEX_ALPHA  = -42.;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  
  if (INDEX_ALPHA != alpha)
  {
    if (table!=0) sm2_free_matrix(table,0, table_N_halo_s-1, 0, table_N_halo_s-1);     
    table = sm2_matrix(0, table_N_halo_s-1, 0, table_N_halo_s-1); 
    double array[3];
    printf("tri_2h_cov %d\n",alpha);
    logsmin = log(limits.halo_s_min);
    logsmax = log(limits.halo_s_max);
    ds = (logsmax - logsmin)/(table_N_halo_s - 1.);
    slog1 = logsmin;
    as = 1./(redshiftshear.zdistrpar_zmax +1.);
    array[2] = alpha*1.0;
    for (i=0; i<table_N_halo_s; i++, slog1+=ds) {  
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<table_N_halo_s; j++, slog2+=ds) {  
	array[1] = exp(slog2);
	//printf("test %le %le\n",l1,l2);
	res =int_GSL_integrate_crude(inner_project_tri_2h_cov,(void*)array,as,0.9999,NULL,1000);
	table[i][j]=log(res);
	table[j][i]=log(res);
//	printf("%le %le %le\n",array[0],array[1],res);

	}
    }
  INDEX_ALPHA=alpha;
  }
  slog1=log(l1);
  slog2=log(l2);
//  printf("tri_2h_cov %le %le\n",l1,l2);
  val = sm2_interpol2d_noextra(table, table_N_halo_s, logsmin, logsmax, ds, slog1, table_N_halo_s, logsmin, logsmax, ds, slog2, 1.0,1.0);
  return exp(val);
}


double project_tri_1h_cov(double l1,double l2,int alpha){
  static double INDEX_ALPHA  = -42.;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  if (INDEX_ALPHA != alpha)
  {
    if (table!=0) sm2_free_matrix(table,0, table_N_halo_s-1, 0, table_N_halo_s-1);     
    table = sm2_matrix(0, table_N_halo_s-1, 0, table_N_halo_s-1); 
    double array[3];
    printf("tri_1h_cov %d\n",alpha);
    logsmin = log(limits.halo_s_min);
    logsmax = log(limits.halo_s_max);
    ds = (logsmax - logsmin)/(table_N_halo_s - 1.);
    slog1 = logsmin;
    as = 1./(redshiftshear.zdistrpar_zmax +1.);
    array[2] = alpha*1.0;
    for (i=0; i<table_N_halo_s; i++, slog1+=ds) {  
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<table_N_halo_s; j++, slog2+=ds) {  
	array[1] = exp(slog2);
	res =int_GSL_integrate_crude(inner_project_tri_1h_cov,(void*)array,as,0.9999,NULL,1000);
	table[i][j]=log(res);
	table[j][i]=log(res);
	//printf("%le %le %le\n",array[0],array[1],res);
	}
    }
  INDEX_ALPHA=alpha;
  }
  slog1=log(l1);
  slog2=log(l2);
  //printf("tri_1h_cov %le %le\n",l1,l2);
  val = sm2_interpol2d_noextra(table, table_N_halo_s, logsmin, logsmax, ds, slog1, table_N_halo_s, logsmin, logsmax, ds, slog2, 1.0,1.0);
  return exp(val);
}


double project_tri_lin_cov(double l1,double l2,int alpha){
  static int INDEX_ALPHA  = -42;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  if (INDEX_ALPHA != alpha)
  {
    if (table!=0) sm2_free_matrix(table,0, table_N_halo_s-1, 0, table_N_halo_s-1);     
    table = sm2_matrix(0, table_N_halo_s-1, 0, table_N_halo_s-1); 
    double array[3];
    //printf("tri_lin_cov %d\n",alpha);

    logsmin = log(limits.halo_s_min);
    logsmax = log(limits.halo_s_max);
    ds = (logsmax - logsmin)/(table_N_halo_s - 1.);
    slog1 = logsmin;
    as = 1./(redshiftshear.zdistrpar_zmax +1.);
    array[2] = alpha*1.0;
    for (i=0; i<table_N_halo_s; i++, slog1+=ds) {  
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<table_N_halo_s; j++, slog2+=ds) {  
	array[1] = exp(slog2);
	//printf("test %le %le\n",l1,l2);
	res =int_GSL_integrate_crude(inner_project_tri_lin_cov,(void*)array,as,0.9999,NULL,1000);
	table[i][j]=log(res);
	table[j][i]=log(res);
	//printf("%le %le %le\n",array[0],array[1],res);
      }
    }
  INDEX_ALPHA=alpha;
  }
  slog1=log(l1);
  slog2=log(l2);
  //printf("tri_lin_cov %.20le %.20le %d\n",l1,l2,alpha);
  val = sm2_interpol2d_noextra(table, table_N_halo_s, logsmin, logsmax, ds, slog1, table_N_halo_s, logsmin, logsmax, ds, slog2,1.0,1.0);
  // printf("tri_lin_cov %le\n",exp(val));
 
  return exp(val);
}






