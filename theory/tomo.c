#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "theory_all.h"
#include "maths.h"
#include "cosmology.h"
#include "EB_functions.h"
#include "tomo.h"

extern con constants;
extern pre precision;
extern lim limits;
extern cosmopara cosmology;
extern configpara configuration;
extern redshiftshearpara redshiftshear;
extern redshiftclusteringpara redshiftclustering;
extern globalpara global;




/*************************************************************************************/
/*      START TOMOGRAPHY ROUTINES                  
 */
/************************************************************************************/


double zdistr_tomo1(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN   = -42.;
  static double ZMAX   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftshear.alpha || BETA_P != redshiftshear.beta_p || Z0 != redshiftshear.z0 || ZMIN !=redshiftshear.zdistrpar_zmin || ZMAX !=redshiftshear.zdistrpar_zmax)
  {
    if(redshiftshear.histogram_zbins != 0 )
    {       //use redshift distribution numeric normalization
    //printf("Using redshift distribution REDSHIFT_FILE %le %le\n",redshiftshear.zdistrpar_zmin,redshiftshear.zdistrpar_zmax);
    norm = sm2_qromb(int_for_zdistr_mock,redshiftshear.zdistrpar_zmin,redshiftshear.zdistrpar_zmax);
    //printf("norm=%le\n",norm);
    }
    
    if((redshiftshear.beta_p>0.) && (redshiftshear.histogram_zbins == 0 ))
    {
      //                      if(redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr,redshiftshear.zdistrpar_zmin,redshiftshear.zdistrpar_zmax));
      }
    }
    ALPHA  = redshiftshear.alpha;
    BETA_P = redshiftshear.beta_p;
    Z0     = redshiftshear.z0;
    ZMIN   = redshiftshear.zdistrpar_zmin;
    ZMAX   = redshiftshear.zdistrpar_zmax;
  }
  
  if(redshiftshear.histogram_zbins != 0)
  {
    if((redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax) && (z>redshiftshear.zdistrpar_zmax || z<redshiftshear.zdistrpar_zmin)) return 0.0;
    //printf("int_for_zdistr_mock %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
  }
  
  if((redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax) && (z>redshiftshear.zdistrpar_zmax || z<redshiftshear.zdistrpar_zmin)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftshear.z0;
  if (fabs(redshiftshear.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftshear.beta_p);
  }
  f=exp(-f);
  return norm*pow(x,redshiftshear.alpha)*f;
}


double zdistr_tomo2(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN2   = -42.;
  static double ZMAX2   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftshear.alpha || BETA_P != redshiftshear.beta_p || Z0 != redshiftshear.z0 || ZMIN2 !=redshiftshear.zdistrpar_zmin2 || ZMAX2 !=redshiftshear.zdistrpar_zmax2)
  {
    if(redshiftshear.histogram_zbins != 0 )
    {       //use redshift distribution numeric normalization
    //printf("Using redshift distribution REDSHIFT_FILE %le %le\n",redshiftshear.zdistrpar_zmin,redshiftshear.zdistrpar_zmax);
    norm = sm2_qromb(int_for_zdistr_mock,redshiftshear.zdistrpar_zmin2,redshiftshear.zdistrpar_zmax2);
    //printf("norm=%le\n",norm);
    }
    
    if((redshiftshear.beta_p>0.) && (redshiftshear.histogram_zbins == 0 ))
    {
      //                      if(redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr,redshiftshear.zdistrpar_zmin2,redshiftshear.zdistrpar_zmax2));
      }
    }
    ALPHA  = redshiftshear.alpha;
    BETA_P = redshiftshear.beta_p;
    Z0     = redshiftshear.z0;
    ZMIN2   = redshiftshear.zdistrpar_zmin2;
    ZMAX2   = redshiftshear.zdistrpar_zmax2;
  }
  if(redshiftshear.histogram_zbins != 0)
  {
    if((redshiftshear.zdistrpar_zmin2 || redshiftshear.zdistrpar_zmax2) && (z>redshiftshear.zdistrpar_zmax2 || z<redshiftshear.zdistrpar_zmin2)) return 0.0;
    //printf("int_for_zdistr_mock %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
  }
  
  if((redshiftshear.zdistrpar_zmin2 || redshiftshear.zdistrpar_zmax2) && (z>redshiftshear.zdistrpar_zmax2 || z<redshiftshear.zdistrpar_zmin2)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftshear.z0;
  if (fabs(redshiftshear.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftshear.beta_p);
  }
  f=exp(-f);
  return norm*pow(x,redshiftshear.alpha)*f;
}




/*************************************************************************************/
/* ! global variable 'global.aglob' is needed here ! */


double int_for_g_tomo1(double aprime)
{
  double chi_glob, chi_prime,val;
  chi_glob = chi(global.aglob);
  chi_prime = chi(aprime);
  val=zdistr_tomo1(1./aprime-1.)*f_K(chi_prime-chi_glob)/f_K(chi_prime)/(aprime*aprime);
  
  return val;
}

double int_for_g_tomo2(double aprime)
{
  double chi_glob, chi_prime,val;
  chi_glob = chi(global.aglob);
  chi_prime = chi(aprime);
  val=zdistr_tomo2(1./aprime-1.)*f_K(chi_prime-chi_glob)/f_K(chi_prime)/(aprime*aprime);
  
  return val;
}

/*redshift weighted lens efficiency factor for a source redshift distribution p_chi(d_chi) WL02(94)*/	
double g_source_tomo(double a, int onetwo)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W_0= -42.;
  static double W_A = -42.;
  
  static double BETA_P  = -42.;
  static double ALPHA  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double Z_MAX2  = -42.;
  static double Z_MIN2   = -42.;
  static double **table;
  
  static double da = 0.0;
  double aa, chi_delta;
  int i;
  /*case of single source redshift at redshiftshear.z0, p(w) = delta(w-w0) */
  if ((redshiftshear.beta_p == 0.0) && (redshiftshear.histogram_zbins == 0)) {
    aa = 1./(1.+redshiftshear.z0);
    if (a<=aa) return 0.0;
    chi_delta = chi(aa);
    return f_K(chi_delta-chi(a))/f_K(chi_delta);
  }
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W_0 != cosmology.w0   || W_A!= cosmology.wa  || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin||Z_MAX2 !=redshiftshear.zdistrpar_zmax2 ||Z_MIN2 !=redshiftshear.zdistrpar_zmin2) 
  {
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_a-1);
    da = (1.-limits.a_min)/(table_N_a-1.);
    table[0][0] = 0.0;
    table[1][0] = 0.0;
    aa = limits.a_min+da;
    /*case of redshift distribution */
    for (i=1;i<table_N_a-1;i++,aa+=da) {
      global.aglob = aa;
      table[0][i] = sm2_qromb(int_for_g_tomo1, limits.a_min, aa);
      table[1][i] = sm2_qromb(int_for_g_tomo2, limits.a_min, aa);
    }
    table[0][table_N_a-1] = 1.;
    table[1][table_N_a-1] = 1.;
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    W_0 = cosmology.w0 ;
    W_A = cosmology.wa ;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   =   redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z_MAX   =   redshiftshear.zdistrpar_zmax;
    Z_MIN   =   redshiftshear.zdistrpar_zmin;
    Z_MAX2   =   redshiftshear.zdistrpar_zmax2;
    Z_MIN2   =   redshiftshear.zdistrpar_zmin2;
  }
  return sm2_interpol(table[onetwo-1], table_N_a, limits.a_min, 1., da, a, 1.0, 1.0);
}


//============================================================	

double int_for_p_shear_shear_tomo(double a)
{
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= g_source_tomo(a,1)/(a*a)*g_source_tomo(a,2)/(a*a)/hoverh0;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; //k in units H0/c
  return res;
}


double P_shear_shear_tomo(double s)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BETA_P  = -42.;
  static double ALPHA  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double Z_MAX2  = -42.;
  static double Z_MIN2  = -42.;
 
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin||Z_MAX2 !=redshiftshear.zdistrpar_zmax2 || Z_MIN2 !=redshiftshear.zdistrpar_zmin2)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_vector(table, 0, table_N_s-1);
    table   = sm2_vector(0, table_N_s-1);
    for (i=0; i<table_N_s; i++, slog+=ds) {
      global.sglob = exp(slog);
      f1 = sm2_qromb(int_for_p_shear_shear_tomo, limits.a_min, 0.7); 
      f2 = sm2_qromo(int_for_p_shear_shear_tomo, .7, 1.0, sm2_midpnt);
      table[i]= log(9./4.*DSQR(cosmology.Omega_m )*2.0*constants.pi_sqr*(f1 + f2)); 
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    Z_MAX2   = redshiftshear.zdistrpar_zmax2;
    Z_MIN2   = redshiftshear.zdistrpar_zmin2;
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


/* ============================================================ *
 * Shear correlation function xi_+ (pm=+1) and xi_- (pm=-1).   *
 * ============================================================ */


double xi_shear_shear_tomo(int pm, double theta)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double SIGMA_8 = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double N_SPEC  = -42.;
  static double ALPHA  = -42.;
  static double BETA_P  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double Z_MAX2  = -42.;
  static double Z_MIN2   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin ||Z_MAX2 !=redshiftshear.zdistrpar_zmax2 ||Z_MIN2 !=redshiftshear.zdistrpar_zmin2 || NONLINEAR != cosmology.nonlinear) 
  {
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    OMB =cosmology.omb;
    H0 =cosmology.h0;
    
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =  redshiftshear.z0;
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    Z_MAX2   = redshiftshear.zdistrpar_zmax2;
    Z_MIN2   = redshiftshear.zdistrpar_zmin2;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_shear_tomo(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}


void xi_via_hankel_shear_shear_tomo(double **xi, double *logthetamin, double *logthetamax)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -42.0, loglmin, dlnl,  lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i, count;
  lP   = fftw_malloc(table_N_thetaH*sizeof(double));
  f_lP = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(table_N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(table_N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-42.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*table_N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = table_N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<table_N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*P_shear_shear_tomo(l);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  for (count=0; count<=1; count++) {
    arg[1] = (count==0 ? 0 : 4);   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<table_N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*table_N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[table_N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<table_N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[count][table_N_thetaH-i-1] = lP[i]/(t*2*constants.pi*table_N_thetaH);
    }
  }
  
  *logthetamin = (nc-table_N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}


//============================================================================
//============================================================================
// SHEAR-MAGNIFICATION ROUTINES
//============================================================================
//============================================================================

double xi_shear_magnification_tomo(int pm, double theta)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double SIGMA_8 = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double N_SPEC  = -42.;
  static double ALPHA  = -42.;
  static double BETA_P  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double Z_MAX2  = -42.;
  static double Z_MIN2   = -42.;
 
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin||Z_MAX2 !=redshiftshear.zdistrpar_zmax2 ||Z_MIN2 !=redshiftshear.zdistrpar_zmin2 || NONLINEAR != cosmology.nonlinear) 
  {
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    SIGMA_8 = cosmology.sigma_8;
    N_SPEC  = cosmology.n_spec;
    OMB =cosmology.omb;
    H0 =cosmology.h0;
    
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =  redshiftshear.z0;
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    Z_MAX2   = redshiftshear.zdistrpar_zmax2;
    Z_MIN2   = redshiftshear.zdistrpar_zmin2;
   
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_magnification_tomo(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return 2.0*res; //as P_shear_shear = 2 P_shear_mag
}


void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -42.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(table_N_thetaH*sizeof(double));
  f_lP = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(table_N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(table_N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-42.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*table_N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = table_N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<table_N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*P_shear_shear_tomo(l);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
   arg[1] = 2;   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<table_N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*table_N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[table_N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<table_N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[0][table_N_thetaH-i-1] = lP[i]/(t*2*constants.pi*table_N_thetaH);
    }
  
  *logthetamin = (nc-table_N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}


//============================================================================
//============================================================================
// SHEAR-POSITION ROUTINES
//============================================================================
//============================================================================

double pf_tomo1(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN   = -42.;
  static double ZMAX   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftclustering.alpha || BETA_P != redshiftclustering.beta_p || Z0 != redshiftclustering.z0 || ZMIN !=redshiftclustering.zdistrpar_zmin || ZMAX !=redshiftclustering.zdistrpar_zmax)
  {
    if(redshiftclustering.histogram_zbins != 0 )
    {       //use redshift distribution numeric normalization
    //printf("Using redshift distribution REDSHIFT_FILE %le %le\n",redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax);
    norm = sm2_qromb(int_for_zdistr_mock,redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax);
    //printf("norm=%le\n",norm);
    }
    
    if((redshiftclustering.beta_p>0.) && (redshiftclustering.histogram_zbins == 0 ))
    {
      //                      if(redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr,redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax));
      }
    }
    ALPHA  = redshiftclustering.alpha;
    BETA_P = redshiftclustering.beta_p;
    Z0     = redshiftclustering.z0;
    ZMIN   = redshiftclustering.zdistrpar_zmin;
    ZMAX   = redshiftclustering.zdistrpar_zmax;
  }
  
  //if single source redshif plane :
  if((redshiftclustering.beta_p==0) && (redshiftclustering.histogram_zbins == 0))
  {
    
    if(fabs(z-redshiftclustering.z0)<1e-2){
      return 1.0;printf("hello\n");}
      else return 0.0;
  }
  
  if(redshiftclustering.histogram_zbins != 0)
  {
    if((redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax) && (z>redshiftclustering.zdistrpar_zmax || z<redshiftclustering.zdistrpar_zmin)) return 0.0;
    //printf("int_for_zdistr_mock %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
  }
  
  if((redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax) && (z>redshiftclustering.zdistrpar_zmax || z<redshiftclustering.zdistrpar_zmin)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftclustering.z0;
  if (fabs(redshiftclustering.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftclustering.beta_p);
  }
  f=exp(-f);
  return norm*pow(x,redshiftclustering.alpha)*f;
}



double pf_tomo2(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN   = -42.;
  static double ZMAX   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftclustering.alpha || BETA_P != redshiftclustering.beta_p || Z0 != redshiftclustering.z0 || ZMIN !=redshiftclustering.zdistrpar_zmin2 || ZMAX !=redshiftclustering.zdistrpar_zmax2)
  {
    if(redshiftclustering.histogram_zbins != 0 )
    {       //use redshift distribution numeric normalization
    //printf("Using redshift distribution REDSHIFT_FILE %le %le\n",redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax);
    norm = sm2_qromb(int_for_zdistr_mock,redshiftclustering.zdistrpar_zmin2,redshiftclustering.zdistrpar_zmax2);
    //printf("norm=%le\n",norm);
    }
    
    if((redshiftclustering.beta_p>0.) && (redshiftclustering.histogram_zbins == 0 ))
    {
      //                      if(redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr,redshiftclustering.zdistrpar_zmin2,redshiftclustering.zdistrpar_zmax2));
      }
    }
    ALPHA  = redshiftclustering.alpha;
    BETA_P = redshiftclustering.beta_p;
    Z0     = redshiftclustering.z0;
    ZMIN   = redshiftclustering.zdistrpar_zmin2;
    ZMAX   = redshiftclustering.zdistrpar_zmax2;
  }
  
  //if single source redshif plane :
  if((redshiftclustering.beta_p==0) && (redshiftclustering.histogram_zbins == 0))
  {
    
    if(fabs(z-redshiftclustering.z0)<1e-2){
      return 1.0;printf("hello\n");}
      else return 0.0;
  }
  
  if(redshiftclustering.histogram_zbins != 0)
  {
    if((redshiftclustering.zdistrpar_zmin2 || redshiftclustering.zdistrpar_zmax2) && (z>redshiftclustering.zdistrpar_zmax2 || z<redshiftclustering.zdistrpar_zmin2)) return 0.0;
    //printf("int_for_zdistr_mock %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
  }
  
  if((redshiftclustering.zdistrpar_zmin2 || redshiftclustering.zdistrpar_zmax2) && (z>redshiftclustering.zdistrpar_zmax2 || z<redshiftclustering.zdistrpar_zmin2)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftclustering.z0;
  if (fabs(redshiftclustering.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftclustering.beta_p);
  }
  f=exp(-f);
  return norm*pow(x,redshiftclustering.alpha)*f;
}



double int_for_p_shear_position_tomo(double a)
{
  double res,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  
  res= g_source_tomo(a,2)*pf_tomo1(1./a-1.)/a/fK;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; 
  return res;
}



double P_shear_position_tomo(double s)  //see Eq. 157 in Schneider 2006 WL
{
  
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BETA_P  = -42.;
  static double ALPHA  = -42.;
  static double Z0      = -42.;
  static double Z1  = -42.;
  static double Z2   = -42.;
  static double BETA_P2  = -42.;
  static double ALPHA2 = -42.;
  static double Z02      = -42.;
  static double Z3 = -42.;
  static double Z4   = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  static double NONLINEAR = -42;
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2,fK,a, k;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z1 !=redshiftshear.zdistrpar_zmin || Z2 !=redshiftshear.zdistrpar_zmin2|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z3 !=redshiftclustering.zdistrpar_zmin || Z4 !=redshiftclustering.zdistrpar_zmin2|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr || NONLINEAR != cosmology.nonlinear)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_vector(table, 0, table_N_s-1);
    table   = sm2_vector(0, table_N_s-1);
  
   for (i=0; i<table_N_s; i++, slog+=ds) {
      global.sglob = exp(slog);
      /*case of single source redshift at redshiftclustering.z0, p(w) = delta(w-w0) */
      if ((redshiftclustering.beta_p == 0.0) && (redshiftclustering.histogram_zbins == 0)) {
	
	a = 1./(1.+redshiftclustering.z0);
	fK     =f_K(chi(a));
	k      = global.sglob/fK;
	table[i]= log(3./2.0*(cosmology.Omega_m)*g_source(a)/a/fK*cosmology.bias*cosmology.rcorr*2.0*constants.pi_sqr*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k); 
	//printf("%le %le %le %le %le %le %le\n",a,g_source(a),hoverh0,fK,Delta_nl_halo(a,k),global.sglob,exp(table[i]));     
      }
      else{
	f1 = sm2_qromb(int_for_p_shear_position_tomo, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_shear_position_tomo, .7, 1.0, sm2_midpnt);
	table[i]= log(3./2.0*(cosmology.Omega_m)*cosmology.bias*cosmology.rcorr*2.0*constants.pi_sqr*(f1 + f2)); 
	//printf("calculating P_shear_position %d\n",i);
      }
    }
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z1   = redshiftshear.zdistrpar_zmin;
    Z2   = redshiftshear.zdistrpar_zmin2;
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z3   = redshiftclustering.zdistrpar_zmin;
    Z4   = redshiftclustering.zdistrpar_zmin2;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


double xi_shear_position_tomo(int pm, double theta)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BETA_P  = -42.;
  static double ALPHA  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double BETA_P2  = -42.;
  static double ALPHA2 = -42.;
  static double Z02      = -42.;
  static double Z_MAX2 = -42.;
  static double Z_MIN2   = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
 {
   OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z_MAX2   = redshiftclustering.zdistrpar_zmax;
    Z_MIN2   = redshiftclustering.zdistrpar_zmin;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;

    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_position(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
    

  return res;
}


void xi_via_hankel_shear_position_tomo(double **xi, double *logthetamin, double *logthetamax)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -42.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(table_N_thetaH*sizeof(double));
  f_lP = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(table_N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(table_N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-42.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*table_N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = table_N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<table_N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*P_shear_position(l);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
    arg[1] = 2;   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<table_N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*table_N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[table_N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<table_N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[0][table_N_thetaH-i-1] = lP[i]/(t*2*constants.pi*table_N_thetaH);
    }
  
  *logthetamin = (nc-table_N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}

double xi_magnification_position_tomo(int pm, double theta)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BETA_P  = -42.;
  static double ALPHA  = -42.;
  static double Z0      = -42.;
  static double Z_MAX  = -42.;
  static double Z_MIN   = -42.;
  static double BETA_P2  = -42.;
  static double ALPHA2 = -42.;
  static double Z02      = -42.;
  static double Z_MAX2 = -42.;
  static double Z_MIN2   = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
 {
   OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   = redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z_MAX2   = redshiftclustering.zdistrpar_zmax;
    Z_MIN2   = redshiftclustering.zdistrpar_zmin;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;

    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_magnification_position(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
    

  return 2.0*res; //as P_shear_position = 2 P_mag_position
}


void xi_via_hankel_magnification_position_tomo(double **xi, double *logthetamin, double *logthetamax)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -42.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(table_N_thetaH*sizeof(double));
  f_lP = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(table_N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(table_N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-42.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*table_N_thetaH-1);
   //dlntheta = (lntmax-lntmin)/(1.0*table_N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = table_N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<table_N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*P_shear_position(l);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
    arg[1] = 0;   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<table_N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*table_N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[table_N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<table_N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[0][table_N_thetaH-i-1] = lP[i]/(t*2*constants.pi*table_N_thetaH);
    }
  
  *logthetamin = (nc-table_N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}

/* ===========================================================================*
 * //Starting routines for position-position power spectrum
 * ==========================================================================*/
double int_for_p_position_position_tomo(double a)
{
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
     
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  
  res= pf(1./a-1.)*pf(1./a-1.)*a*a*hoverh0/fK/fK;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; 
  return res;
}



double P_position_position_tomo(double s)  //see Eq. 157 in Schneider 2006 WL
{
  
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  static double NONLINEAR = -42;
  
  static double BETA_P2  = -42.;
  static double ALPHA2 = -42.;
  static double Z02      = -42.;
  static double Z_MAX2 = -42.;
  static double Z_MIN2   = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2,fK,a, hoverh0,k;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 ||  BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_vector(table, 0, table_N_s-1);
    table   = sm2_vector(0, table_N_s-1);
  
    for (i=0; i<table_N_s; i++, slog+=ds) {
      global.sglob = exp(slog);
      /*case of single source redshift at redshiftclustering.z0, p(w) = delta(w-w0) */
      if ((redshiftclustering.beta_p == 0.0) && (redshiftclustering.histogram_zbins == 0)) {
	
	a = 1./(1.+redshiftclustering.z0);
	fK     =f_K(chi(a));
	k      = global.sglob/fK;
	hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
	table[i]= log(cosmology.bias*cosmology.bias*a*a*hoverh0/fK/fK*2.0*constants.pi_sqr*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k); 
	//printf("%le %le %le %le %le  %le\n",a,hoverh0,fK,Delta_nl_halo(a,k),global.sglob,exp(table[i]));
	
      }
      else{
	f1 = sm2_qromb(int_for_p_position_position_tomo, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_position_position_tomo, .7, 1.0, sm2_midpnt);
	table[i]=log(cosmology.bias*cosmology.bias*2.0*constants.pi_sqr*(f1 + f2)); 
	//printf("calculating P_shear_position %d\n",i);
      }
    } 
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z_MAX2   = redshiftclustering.zdistrpar_zmax;
    Z_MIN2   = redshiftclustering.zdistrpar_zmin;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}



double xi_position_position_tomo(int pm, double theta)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  
  static double BETA_P2  = -42.;
  static double ALPHA2 = -42.;
  static double Z02      = -42.;
  static double Z_MAX2 = -42.;
  static double Z_MIN2   = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 ||  BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
  {
    OMEGA_M =   cosmology.Omega_m ;
    OMEGA_V =   cosmology.Omega_v ;
    N_SPEC  =   cosmology.n_spec;
    OMB = cosmology.omb;
    H0 = cosmology.h0;
    SIGMA_8 =  cosmology.sigma_8;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    NONLINEAR = cosmology.nonlinear;
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z_MAX2   = redshiftclustering.zdistrpar_zmax;
    Z_MIN2   = redshiftclustering.zdistrpar_zmin;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_position_position_tomo(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}



void xi_via_hankel_position_position_tomo(double **xi, double *logthetamin, double *logthetamax)
{
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -42.0, loglmin, dlnl, lnrc, arg[2];
  static int nc;
  
  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(table_N_thetaH*sizeof(double));
  f_lP = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((table_N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(table_N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(table_N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-42.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*table_N_thetaH-1);
//     lntmin   = log(limits.xi_via_hankel_theta_min);
//     lntmax   = log(limits.xi_via_hankel_theta_max);
//     dlntheta = (lntmax-lntmin)/(1.0*table_N_thetaH-1);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = table_N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<table_N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*P_position_position_tomo(l);
  }
  
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
    arg[1] = 0;   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<table_N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*table_N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[table_N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<table_N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[0][table_N_thetaH-i-1] = lP[i]/(t*2*constants.pi*table_N_thetaH);
    }
  
  
  *logthetamin = (nc-table_N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}


