#include <math.h>
#include <stdlib.h>
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
#include "halo.h"

extern con constants;
extern pre precision;
extern lim limits;
extern cosmopara cosmology;
extern configpara configuration;
extern redshiftpara redshift;
extern globalpara global;
extern tomopara tomo;



/*************************************************************************************/
/*      START TOMOGRAPHY ROUTINES                  
 */
/************************************************************************************/

double indexcalc(int a,int b)
{
  int i,j,k;
  static int NEW=-42;
  static double **table;
  
  if (NEW!=0){
    if (table!=0) sm2_free_matrix(table, 0, tomo.shear_Nbin-1, 0, tomo.shear_Nbin-1);
    table = sm2_matrix(0, tomo.shear_Nbin-1, 0,tomo.shear_Nbin-1);
    k=0;
    for (i=0; i<tomo.shear_Nbin;i++){
     for (j=i; j<tomo.shear_Nbin;j++,k++){
      table[i][j]=k;
      table[j][i]=k;       
      }
    }
    NEW=0;
  } 
  return table[a][b];
}

double zdistr_tomo(double z, int i)
{
    double static norm;
    int static normi=-42;
  //printf("%le %le\n",tomo.shear_zmax[i],tomo.shear_zmin[i]);
  //First, compute the normalization
  if(normi !=i){  
  norm = sm2_qromb(int_for_zdistr_mock,tomo.shear_zmin[i],tomo.shear_zmax[i]);
   // printf("norm=%le\n",norm);
  normi=i;
  }
    if(z>tomo.shear_zmax[i] || z<tomo.shear_zmin[i]) return 0.0;
    //printf("zdistr_tomo %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
}

double zdistr_tomo2(double z, int i)
{
    double static norm;
    int static normi=-42;
  //printf("%le %le\n",tomo.shear_zmax[i],tomo.shear_zmin[i]);
  //First, compute the normalization
  if(normi !=i){  
  norm = sm2_qromb(int_for_zdistr_mock,tomo.shear_zmin[i],tomo.shear_zmax[i]);
    //printf("norm=%le\n",norm);
  normi=i;
  }
    if(z>tomo.shear_zmax[i] || z<tomo.shear_zmin[i]) return 0.0;
    //printf("zdistr_tomo2 %le %le %le %le\n",z,int_for_zdistr_mock(z),norm,int_for_zdistr_mock(z)/norm);
    return int_for_zdistr_mock(z)/norm;
}


/*************************************************************************************/
/* ! global variable 'global.aglob' is needed here ! */


double int_for_g_tomo(double aprime,void *params)
{
  double chi_glob, chi_prime,val;
  double *ar = (double *) params;
  int zbin= (int) ar[0];
  //if (global.aglob < 1./(tomo.shear_zmax[zbin]+1.) ) return 0.0;
      
  chi_glob = chi(global.aglob);
  chi_prime = chi(aprime);
  val=zdistr_tomo(1./aprime-1.,zbin)*f_K(chi_prime-chi_glob)/f_K(chi_prime)/(aprime*aprime);
 // printf("tomo %d  %le %le\n",zbin,1./aprime-1., zdistr_tomo(1./aprime-1.,(int) ar[0]));
  return val;
}

/*redshift weighted lens efficiency factor for a source redshift distribution p_chi(d_chi) WL02(94)*/	
double g_source_tomo(double a, int zbin)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W_0= -42.;
  static double W_A = -42.;
  
  static double **table;
  
  static double da = 0.0;
  double aa;
  int i,j;
  double array[1];
 /*case of single source redshift at redshift.shear_z0, p(w) = delta(w-w0) */
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W_0 != cosmology.w0   || W_A!= cosmology.wa) 
  {
    if (table!=0) sm2_free_matrix(table, 0, tomo.shear_Nbin-1, 0, table_N_a-1);
    table   = sm2_matrix(0, tomo.shear_Nbin-1, 0, table_N_a-1);
    da = (0.999999-1./(redshift.shear_zdistrpar_zmax+1.))/(table_N_a-1);
    
    for (j=0;j<tomo.shear_Nbin;j++) {
      array[0]=(double) j;
      aa = 1./(redshift.shear_zdistrpar_zmax+1.);
    /*case of redshift distribution */
    for (i=0;i<table_N_a;i++,aa+=da) {
      global.aglob = aa;
      table[j][i] = int_GSL_integrate_qag2(int_for_g_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),aa,NULL,4000);
      //printf(" tomo  %le %le %le\n",aa,1./aa-1,table[j][i]);
     }
//     if (j==4){
// 	FILE *F;
// 	F=fopen("../../../data/baryons/g_source_tomo","w");
// 	aa = 1./(redshift.shear_zdistrpar_zmax+1.);
// 	for (i=0;i<table_N_a;i++,aa+=da) {
// 	fprintf(F,"%le %le %le %le %le %le\n",1./aa-1,table[0][i],table[1][i],table[2][i],table[3][i],table[4][i]);
//       }
//       fclose(F);
//     }
    
    }

    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    W_0 = cosmology.w0 ;
    W_A = cosmology.wa ;
  }
  if (a<1./(redshift.shear_zdistrpar_zmax+1.)) return 0.0;
  return sm2_interpol(table[zbin], table_N_a, 1./(redshift.shear_zdistrpar_zmax+1.), 0.999999, da, a, 1.0, 1.0);
}


//============================================================	

double int_for_p_shear_shear_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= g_source_tomo(a,(int) ar[0])/(a*a)*g_source_tomo(a,(int) ar[1])/(a*a)/hoverh0;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; //k in units H0/c
  return res;
}


double P_shear_shear_tomo(double s, int Npowerspec)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
 
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog;
  int i,j,k,l;
  
  double array[2],res = 0.0;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
   
    if (table!=0) sm2_free_matrix(table, 0, tomo.shear_Npowerspectra-1, 0, table_N_s-1);
    table   = sm2_matrix(0, tomo.shear_Npowerspectra-1, 0, table_N_s-1);
   
    l=0;
    for (k=0; k<tomo.shear_Nbin; k++) {
      array[0]=(double) k;
      for (j=k; j<tomo.shear_Nbin; j++,l++) {
	array[1]=(double) j;
	 slog = logsmin;
	 for (i=0; i<table_N_s; i++, slog+=ds) {
	  global.sglob = exp(slog);
	  res =int_GSL_integrate_crude(int_for_p_shear_shear_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),0.999999,NULL,1000);
	  table[l][i]= log(9./4.*DSQR(cosmology.Omega_m )*2.0*constants.pi_sqr*(res)); 
	  //printf("%d %d %d %le %le\n",k,j,l,exp(slog),table[l][i]);
	}
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
  }
  slog = log(s);
  f1 = sm2_interpol(table[Npowerspec], table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


/* ============================================================ *
 * Shear correlation function xi_+ (pm=+1) and xi_- (pm=-1).   *
 * ============================================================ */


double xi_shear_shear_plus_tomo(int pm, double theta,int Npowerspec)
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
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || BETA_P != redshift.shear_beta_p || ALPHA != redshift.shear_alpha || Z0 != redshift.shear_z0 || NONLINEAR != cosmology.nonlinear) 
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
    
    BETA_P  =   redshift.shear_beta_p;
    ALPHA   = redshift.shear_alpha;
    Z0      =  redshift.shear_z0;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_shear_tomo(table, &logthetamin, &logthetamax,Npowerspec);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}


void xi_via_hankel_shear_shear_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec)
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
    lP[i] = l*P_shear_shear_tomo(l,Npowerspec);
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

double xi_shear_magnification_tomo(int pm, double theta,int Npowerspec)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double SIGMA_8 = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double N_SPEC  = -42.;  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || NONLINEAR != cosmology.nonlinear) 
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
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_magnification_tomo(table, &logthetamin, &logthetamax,Npowerspec);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return 2.0*res; //as P_shear_shear = 2 P_shear_mag
}


void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec)
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
    lP[i] = l*P_shear_shear_tomo(l,Npowerspec);
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


double int_for_p_shear_position_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  
  res= g_source_tomo(a,(int) ar[1])*zdistr_tomo(1./a-1.,(int) (ar[0]))/a/fK; //ar[0]<ar[1] otherwise no sp-correlation; foreground position is correlated with background sources but not vv 
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; 
 return res;
}



double P_shear_position_tomo(double s, int Npowerspec)  //see Eq. 157 in Schneider 2006 WL
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  static double NONLINEAR = -42;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double slog, k,f1;
  int i,j,l;
  double array[2],res = 0.0;
   
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS != cosmology.bias ||RCORR !=cosmology.rcorr || NONLINEAR != cosmology.nonlinear)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_matrix(table, 0, tomo.shear_Npowerspectra-1, 0, table_N_s-1);
    table   = sm2_matrix(0, tomo.shear_Npowerspectra-1, 0, table_N_s-1);
  
    l=0;
    for (k=0; k<tomo.shear_Nbin; k++) {
      array[0]=(double) k;
      for (j=k; j<tomo.shear_Nbin; j++,l++) {
	array[1]=(double) j;
	 slog = logsmin;
	 for (i=0; i<table_N_s; i++, slog+=ds) {
	  global.sglob = exp(slog);
	  res =int_GSL_integrate_crude(int_for_p_shear_position_tomo,(void*)array,limits.a_min,0.9999,NULL,1000);
	  table[l][i]= log(3./2.0*(cosmology.Omega_m)*cosmology.bias*cosmology.rcorr*2.0*constants.pi_sqr*res); 
	  //printf("%d %d %d %le %le\n",k,j,l,exp(slog),table[l][i]);
	}
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
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
  }
  slog = log(s);
  f1 = sm2_interpol(table[Npowerspec], table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


double xi_shear_position_tomo(int pm, double theta,int Npowerspec)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
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
    
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_position_tomo(table, &logthetamin, &logthetamax,Npowerspec);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  
  
  return res;
}


void xi_via_hankel_shear_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec)
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
    lP[i] = l*P_shear_position_tomo(l,Npowerspec);
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

double xi_magnification_position_tomo(int pm, double theta,int Npowerspec)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
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
    
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_magnification_position_tomo(table, &logthetamin, &logthetamax,Npowerspec);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  
  
  return 2.0*res; //as P_shear_position = 2 P_mag_position
}


void xi_via_hankel_magnification_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec)
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
    lP[i] = l*P_shear_position_tomo(l,Npowerspec);
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
double int_for_p_position_position_tomo(double a, void *params)
{
  double res,hoverh0,ell, fK, k;
  double *ar = (double *) params;
 if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  
  res= zdistr_tomo(1./a-1.,(int) ar[0])*zdistr_tomo2(1./a-1.,(int) ar[1])*a*a*hoverh0/fK/fK;
  //printf("%le %d %d %le %le\n",1./a-1.,(int) ar[0],(int) ar[1],zdistr_tomo(1./a-1.,(int) ar[0]),zdistr_tomo2(1./a-1.,(int) ar[1]));
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; 
 // printf("res integrand %le %le\n",a,Delta_nl_halo(a,k/cosmology.coverH0));
  return res;
}



double P_position_position_tomo(double s, int Npowerspec)  //see Eq. 157 in Schneider 2006 WL
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
  
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,k;
  int i,j,l;
  
  double array[2],res = 0.0;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    
    if (table!=0) sm2_free_matrix(table, 0, tomo.shear_Npowerspectra-1, 0, table_N_s-1); 
    table   = sm2_matrix(0, tomo.shear_Npowerspectra-1, 0, table_N_s-1);
    
    l=0;
    for (k=0; k<tomo.shear_Nbin; k++) {
      array[0]=(double) k;
      for (j=k; j<tomo.shear_Nbin; j++,l++) {
	array[1]=(double) j;
	 slog = logsmin;
	 for (i=0; i<table_N_s; i++, slog+=ds) {
	  global.sglob = exp(slog);
	  res =int_GSL_integrate_crude(int_for_p_position_position_tomo,(void*)array,limits.a_min,0.9999,NULL,1000);
	  table[l][i]= log(cosmology.bias*cosmology.bias*2.0*constants.pi_sqr*res); 
	  //printf("%d %d %d %le %le\n",k,j,l,exp(slog),table[l][i]);
	}
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
    
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
  }
  slog = log(s);
  f1 = sm2_interpol(table[Npowerspec], table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}



double xi_position_position_tomo(int pm, double theta,int Npowerspec)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0= -42.;
  static double WA = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0   = -42.;
  static double SIGMA_8 = -42.;
  
  static double BIAS   = -42.;
  static double RCORR   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear)
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
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_position_position_tomo(table, &logthetamin, &logthetamax,Npowerspec);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}



void xi_via_hankel_position_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec)
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
    lP[i] = l*P_position_position_tomo(l,Npowerspec);
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


