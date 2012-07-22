#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include "theory_all.h"
#include "maths.h"
#include "cosmology.h"
#include "fftw3.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

//====
//in 1/Mpc 
#define k_max_coyote 2.416000 
#define k_min_coyote 0.004833 
//#define k_min_coyote 0.008371 0.006835
//===
#define a_max_coyote 1.0
#define a_min_coyote 0.5

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


extern con constants;
extern pre precision;
extern lim limits;
extern cosmopara cosmology;
extern configpara configuration;
extern redshiftshearpara redshiftshear;
extern redshiftclusteringpara redshiftclustering;
extern globalpara global;
extern nuisancepara nuisance;


gsl_interp_accel *zaccel5;
gsl_spline *redshift_distrib_spline;

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//c evolution of omega matter and omega lamda with expansion factor

void omega_a(double aa,double *om_m,double *om_v)
{  
  double a2,omega_curv;
  a2=aa*aa;
  omega_curv=1.0-cosmology.Omega_m- cosmology.Omega_v;
  *om_m=cosmology.Omega_m /(cosmology.Omega_m +aa*(omv_vareos(aa) *a2 +omega_curv));
  *om_v=omv_vareos(aa)*a2*aa/(cosmology.Omega_m+aa*(a2*omv_vareos(aa) +omega_curv)); 
}

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Dark energy routines
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//variable Omega_v
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double omv_vareos(double a)
{
  return(cosmology.Omega_v*exp(-3.*((cosmology.w0+cosmology.wa+1.)*log(a)+cosmology.wa*(1.-a))));
}

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//growth factor including Dark energy parameters w0, wa
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double growfac(double a, double omm, double omv)
{
  //variables omm,omv ignored!
  const double MINA=1.e-8;
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double ai[table_N_a];
  static double table[table_N_a];
  double res;
  
  gsl_interp *intf=gsl_interp_alloc(gsl_interp_linear,table_N_a);
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  
  if (OMEGA_M != cosmology.Omega_m || OMEGA_V != cosmology.Omega_v || W0 != 
    cosmology.w0 || WA != cosmology.wa)
  {
    OMEGA_M = cosmology.Omega_m;
    OMEGA_V = cosmology.Omega_v;
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    
    int i;
    const gsl_odeiv_step_type *T=gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,2);
    gsl_odeiv_control *c=gsl_odeiv_control_y_new(1.e-6,0.0);  
    //abs./rel. allowed error on y
    gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(2);
    
    double t=MINA;            //start a
    double t1=1.;                //final a
    double h=1.e-6;              //initial step size
    double y[2]={MINA,MINA};   //initial conditions
    double norm;
    double par[4]={OMEGA_M,OMEGA_V,W0,WA};    
    //Omega_m,Omega_v,w0,wa
    gsl_odeiv_system sys={func_for_growfac,NULL,2,&par};
    
    for (i=1;i<=table_N_a;i++) {
      ai[i-1]=i*t1/(1.*table_N_a);
      while(t<ai[i-1]) 
	gsl_odeiv_evolve_apply(e,c,s,&sys,&t,ai[i-1],&h,y);
      if (i==1) norm=y[0]/ai[i-1];
      table[i-1]=y[0]/norm;
    }
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
  }
  gsl_interp_init(intf,ai,table,table_N_a);
  res=gsl_interp_eval(intf,ai,table,a,acc);
  gsl_interp_accel_free(acc);
  gsl_interp_free(intf);
  return(res);
}

//function for growfac (DGL)
int func_for_growfac(double a,const double y[],double f[],void *params)
{
  double *p=(double *)params;
  if (a == 0) {
    printf("a=0 in function 'func_for_growfac'!\n");
    exit(1);
  }
  double aa=a*a;
  double omegam=p[0]/(aa*a);
  double omegav=p[1]*exp(-3.*((p[2]+p[3]+1)*log(a)+p[3]*(1.-a)));
  double hub=omegam+(1-p[0]-p[1])/(a*a)+omegav;
  f[0]=y[1];
  
  f[1]=y[0]*3.*p[0]/(2.*hub*aa*aa*a)-y[1]/a*(2.-(omegam+(3.*(p[2]+p[3]*(1.-a))+1)*omegav)/(2.*hub));
  return GSL_SUCCESS;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// BEGIN EH wiggled routines
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/* ---------------------------- TFmdm_onek_hmpc() ---------------------- */

double Tsqr_EH_wiggle(double kk)
/* Given a wavenumber in h Mpc^-1, return the transfer function for the
 *cosmology held in the global variables. */
/* Input: kk -- Wavenumber in h Mpc^-1 */
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double OMB   = -42.;
  static double H0  = -42.;
  
  /* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
  /* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
  double omega_hdm=0.0; // will be changed later
  int degen_hdm=1;
  
  
  /* The following are set in TFmdm_set_cosm() */
  double  static alpha_gamma,	/* sqrt(alpha_nu) */
  alpha_nu,	/* The small-scale suppression */
  beta_c,		/* The correction to the log in the small-scale */
  num_degen_hdm,	/* Number of degenerate massive neutrino species */
  f_baryon,	/* Baryon fraction */
  f_bnu,		/* Baryon + Massive Neutrino fraction */
  f_cb,		/* Baryon + CDM fraction */
  f_cdm,		/* CDM fraction */
  f_hdm,		/* Massive Neutrino fraction */
  growth_k0,	/* D_1(z) -- the growth function as k->0 */
  growth_to_z0,	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
  k_equality,	/* The comoving wave number of the horizon at equality*/
  obhh,		/* Omega_baryon * hubble^2 */
  omega_curv,	/* = 1 - omega_matter - omega_lambda */
  omega_lambda_z, /* Omega_lambda at the given redshift */
  omega_matter_z,	/* Omega_matter at the given redshift */
  omhh,		/* Omega_matter * hubble^2 */
  onhh,		/* Omega_hdm * hubble^2 */
  p_c,		/* The correction to the exponent before drag epoch */
  p_cb,		/* The correction to the exponent after drag epoch */
  sound_horizon_fit,  /* The sound horizon at the drag epoch */
  theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
  y_drag,		/* Ratio of z_equality to z_drag */
  z_drag,		/* Redshift of the drag epoch */
  z_equality;	/* Redshift of matter-radiation equality */
  
  /* The following are set in TFmdm_onek_mpc() */
  double static gamma_eff,	/* Effective \Gamma */
  growth_cb,	/* Growth factor for CDM+Baryon perturbations */
  growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
  max_fs_correction,  /* Correction near maximal free streaming */
  qq,		/* Wavenumber rescaled by \Gamma */
  qq_eff,		/* Wavenumber rescaled by effective Gamma */
  qq_nu,		/* Wavenumber compared to maximal free streaming */
  tf_master,	/* Master TF */
  tf_sup,		/* Suppressed TF */
  y_freestream; 	/* The epoch of free-streaming for a given scale */
  
  /* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
  double  static tf_cb,		/* The transfer function for density-weighted
  CDM + Baryon perturbations. */
  tf_cbnu;	/* The transfer function for density-weighted
  CDM + Baryon + Massive Neutrino perturbations. */
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || OMB !=cosmology.omb || H0 !=cosmology.h0)
  {
    OMEGA_M = cosmology.Omega_m;
    OMEGA_V =cosmology.Omega_v;
    OMB=cosmology.omb;
    H0= cosmology.h0;
  
    
    /* -------------------------Beginning of Wayne Hu EH99 routine ------------------------ */
    /* -------------------------- TFmdm_set_cosm() ------------------------------------------------- */
    
    double z_drag_b1, z_drag_b2, omega_denom;
    
    theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */
    
    if (degen_hdm<1) degen_hdm=1;
    num_degen_hdm = (double) degen_hdm;	
    /* Have to save this for TFmdm_onek_mpc() */
    /* This routine would crash if baryons or neutrinos were zero, 
     *	so don't allow that */
    if (cosmology.omb<=0) cosmology.omb=1e-5;
    if (omega_hdm<=0) omega_hdm=1e-5;
    
    omega_curv = 1.0-cosmology.Omega_m-cosmology.Omega_v;
    omhh = cosmology.Omega_m*SQR(cosmology.h0);
    obhh = cosmology.omb*SQR(cosmology.h0);
    onhh = omega_hdm*SQR(cosmology.h0);
    f_baryon = cosmology.omb/cosmology.Omega_m;
    f_hdm = omega_hdm/cosmology.Omega_m;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;
    
    /* Compute the equality scale. */
    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
    k_equality = 0.0746*omhh/SQR(theta_cmb);
    
 //!!!!!!!! added this line as Transfer function is calculated at time of matter-rad equality
double  redshift=z_equality;  
 //!!!!!!!!
    /* Compute the drag epoch and sound horizon */
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2));
    y_drag = z_equality/(1.0+z_drag);
    
    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));
    
    /* Set up for the free-streaming & infall growth function */
    p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));
    
    omega_denom = cosmology.Omega_v+SQR(1.0+redshift)*(omega_curv+cosmology.Omega_m*(1.0+redshift));
    omega_lambda_z = cosmology.Omega_v/omega_denom;
    omega_matter_z = cosmology.Omega_m*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
    growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)-omega_lambda_z+ (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
    growth_to_z0 = z_equality*2.5*cosmology.Omega_m/(pow(cosmology.Omega_m,4.0/7.0) -cosmology.Omega_v + (1.0+cosmology.Omega_m/2.0)*(1.0+cosmology.Omega_v/70.0));
    growth_to_z0 = growth_k0/growth_to_z0;	
    
    /* Compute small-scale suppression */
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*pow(1+y_drag,p_cb-p_c)*(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    /* Done setting scalar variables */
  }
  
  kk=kk*cosmology.h0;
  
  double tf_sup_L, tf_sup_C;
  double temp1, temp2;
  
  qq = kk/omhh*SQR(theta_cmb);
  
  /* Compute the scale-dependent growth functions */
  y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*SQR(num_degen_hdm*qq/f_hdm);
  temp1 = pow(growth_k0, 1.0-p_cb);
  temp2 = pow(growth_k0/(1+y_freestream),0.7);
  growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
  growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;
  
  /* Compute the master function */
  gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/(1+SQR(SQR(kk*sound_horizon_fit*0.43))));
  qq_eff = qq*omhh/gamma_eff;
  
  tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
  tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
  tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));
  
  qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
  max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/	(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
  tf_master = tf_sup*max_fs_correction;
  
  /* Now compute the CDM+HDM+baryon transfer functions */
  tf_cb = tf_master*growth_cb/growth_k0;
  tf_cbnu = tf_master*growth_cbnu/growth_k0;
  return tf_cb*tf_cb*pow(kk,cosmology.n_spec);
}


double int_for_sigma_8_EH_wiggle(double k, void *args)
{
  double 	kR, res, x;
  kR = k/375.;
  x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);
  res = DSQR(k)*Tsqr_EH_wiggle(k/cosmology.coverH0)*x*x;
  return res;
}


double sigma_8_sqr_EH_wiggle()    
{
  static double N_SPEC = -42.;
  static double res = -42.;
  static double OMEGA_M = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  double integral;
  
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| N_SPEC != cosmology.n_spec)
  {
    // full integral - method of choice!
    integral = int_GSL_integrate_qag(int_for_sigma_8_EH_wiggle,NULL,1e-4,1e6,NULL,2048);
    res = 4.5/constants.pi/constants.pi*integral;   //see PD97, eq. 29
    
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    N_SPEC = cosmology.n_spec;
  }
  assert(res>0.0);
  return res;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// END EH wiggled routines
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//c Efstathiou, Bond, White (1992) transfer function
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//linear scale-free power spectrum 
double Delta_L(double k)
{
  
  if (configuration.transferfunction_EH99==1){
    double kkorr=k*cosmology.coverH0;
//     printf("Entering\n");
    double k3=kkorr*kkorr*kkorr;
    double fac=cosmology.sigma_8*cosmology.sigma_8*k3/(2.*constants.pi*constants.pi*sigma_8_sqr());
   // printf("entering eisenstein transfer function\n");
    return(fac*Tsqr_EH(kkorr));
  }    //EH 99
  if (configuration.transferfunction_EH99==2){
    double kkorr=k*cosmology.coverH0;
   // printf("KKORR=%le %le %le\n",kkorr,Tsqr_EH_wiggle(kkorr/3000.),Tsqr_EH(kkorr));
    double k3=kkorr*kkorr*kkorr;
    double  fac=cosmology.sigma_8*cosmology.sigma_8*k3/(2.*constants.pi*constants.pi*sigma_8_sqr_EH_wiggle());
    //printf("entering wiggled transfer function\n");
    return(fac*Tsqr_EH_wiggle(kkorr/cosmology.coverH0));
  }    //EH 99
  else {      //original EBW 92
    double keff,q,q8,tk,tk8;
    //printf("entering efstathiou  transfer function\n");
    double Gamma=cosmology.Omega_m*cosmology.h0*exp(-cosmology.omb*(1.+sqrt(2.0*cosmology.h0)/cosmology.Omega_m));
    keff=0.172+0.011*log(Gamma/0.36)*log(Gamma/0.36);
    q=1e-20+k/Gamma;
    q8=1e-20+keff/Gamma;
    tk=1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+(1.7*q)*(1.7*q),1.13),(1/1.13));
    tk8=1/pow(1+pow(6.4*q8+pow(3.0*q8,1.5)+(1.7*q8)*(1.7*q8),1.13),(1/1.13));
    
    return cosmology.sigma_8*cosmology.sigma_8*pow(q/q8,3+cosmology.n_spec)*tk*tk/tk8/tk8;
  }
}


double Delta_L_tab(double rk)
{
  static double OMEGA_M = -42.;
  static double N_SPEC  = -42.;
  static double SIGMA_8 = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  
  static double table_P[table_N_klin];
  static double dk = .0, logkmin = .0, logkmax = .0;
  
  double klog,f1;
  int i;     
  
  if (rk < limits.Delta_L_klin_min || rk > limits.Delta_L_klin_max){ 
    return Delta_L(rk);
    //printf("outside Delta_L_tab\n");   
  }
  else{  
    if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec)
    {
      SIGMA_8 = cosmology.sigma_8;
      N_SPEC  = cosmology.n_spec;
      H0=cosmology.h0;
      OMB=cosmology.omb;
      OMEGA_M=cosmology.Omega_m;
      
      logkmin = log(limits.Delta_L_klin_min);
      logkmax = log(limits.Delta_L_klin_max);
      dk = (logkmax - logkmin)/(table_N_klin-1.);
      klog = logkmin;
      
      for (i=0; i<table_N_klin; i++, klog += dk) {  
	table_P[i]=log(Delta_L(exp(klog)));
      }
      // printf("finished Delta_L_tab\n");   
    }
  }
  klog=log(rk);
  f1=sm2_interpol(table_P, table_N_klin, logkmin, logkmax, dk,klog, 1.0,1.0 );
  return exp(f1);  
}


double p_lin(double k_NL,double a){     
     static double OMEGA_M = -42.;
     static double OMEGA_V = -42.;
     static double N_SPEC  = -42.;
     static double OMB   = -42.;
     static double H0  = -42.;
     static double SIGMA_8 = -42.;
     static double W0      = -42.;
     static double WA      = -42.;
     static double **table_P_Lz = 0;
     static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
       
     double plin,rk,om_m,om_v,amp,ampsqr,grow0,aa,klog,Pdelta,val;
     
//     printf("plin test a=%le k=%le\n",a,k_NL);
     int i,j;
     if (a >= 0.99999){a =0.99999;}
     if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  ||W0 != cosmology.w0 || WA != cosmology.wa || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || SIGMA_8 != cosmology.sigma_8) {
     if (table_P_Lz!=0) sm2_free_matrix(table_P_Lz,0, table_N_a-1, 0, table_N_k-1);  
     table_P_Lz = sm2_matrix(0, table_N_a-1, 0, table_N_k-1);
     //printf("new P_Delta\n");
     grow0=growfac(1.,1.,1.);
     da = (1. - limits.a_min)/(table_N_a-1.);
     aa = limits.a_min;
     for (i=0; i<table_N_a; i++, aa +=da) {
      if(aa>1.0) aa=1.0;
      omega_a(aa,&om_m,&om_v);
      amp=growfac(aa,1.,1.)/grow0;
      ampsqr=amp*amp;
      
      logkmin = log(limits.Pdelta_halo_k_min);
      logkmax = log(limits.Pdelta_halo_k_max);
      dk = (logkmax - logkmin)/(table_N_k-1.);
      klog = logkmin;
        for (j=0; j<table_N_k; j++, klog += dk) {
	 rk=exp(klog)/cosmology.coverH0;
     //    printf("plin test a=%le k=%le\n",aa,rk);
	 plin=ampsqr*Delta_L_tab(rk);
//	 printf("plin test a=%le k=%le\n",aa,rk); 
	  Pdelta=plin;
	//table_k[j] = klog;
	table_P_Lz[i][j] = log(2*constants.pi_sqr*Pdelta) - 3.*klog;
	//printf("%le %le\n",rk,table_P_Lz[i][j]);
	}
	//printf("%le %le\n",limits.Pdelta_halo_k_min,limits.Pdelta_halo_k_max);
//	printf("table finished\n");
//	upper = cosmology.n_spec-4.0;
        //sm2_spline(table_k-1, table_P-1, table_N_k, cosmology.n_spec, 0.0, y2-1);
// 	klog = logkmin;
// 	for (j=0; j<table_N_k; j++, klog += dk) {
// 	    sm2_splint(table_k-1, table_P-1, y2-1, table_N_k, klog, &val);
// 	    table_P_Lz[i][j] = val;
// 	    }
      }
//      printf("table finished\n");
      OMEGA_M = cosmology.Omega_m ;
      OMEGA_V = cosmology.Omega_v ;
      N_SPEC  = cosmology.n_spec;
      OMB = cosmology.omb;
      H0=cosmology.h0; 
      SIGMA_8 = cosmology.sigma_8; 
      W0      = cosmology.w0;
      WA      = cosmology.wa;
      //printf("P_delta_table finished\n");
   }    
   klog = log(k_NL); 
   //printf("%le\n",exp(klog));
//   upper = cosmology.n_spec-4.0;
   
   val = sm2_interpol2d(table_P_Lz, table_N_a, limits.a_min, 1., da, a, table_N_k, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
   //if(k_NL>1.e3) printf("%le %le\n",k_NL,exp(val));
   return exp(val);     
}



//c /*Eisenstein &Hu no-wiggle transfer function*/  /* from Martin White */
double Tsqr_EH(double k)
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;
  
  /* other input parameters */
  hubble = cosmology.h0;
  
  omegam = cosmology.Omega_m;
  ombh2 = cosmology.omb* cosmology.h0 *cosmology.h0 ;
  
  if(cosmology.omb == 0)
    ombh2 = 0.04 * cosmology.h0 *cosmology.h0 ;
  
  k *= 1./cosmology.coverH0;//(3.085678e24 / UnitLength_in_cm);       /* convert to h/Mpc */
  
  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * 
  log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
  + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * 
  (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma; 
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp*tmp)*pow(k,cosmology.n_spec);
}

//Further routines needed for Eisenstein transfer function
double int_for_sigma_8(double k, void *args)
{
  double 	kR, res, x;
  kR = k/375.;
  x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);
  res = DSQR(k)*Tsqr_EH(k)*x*x;
  return res;
}

double int_GSL_integrate_qag(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
  gsl_integration_qag(&F,a,b,1e-5,1e-5,niter,GSL_INTEG_GAUSS41,w,&res,&err);
  if(NULL!=error)
    *error=err;
  gsl_integration_workspace_free(w);
  return res;
}


double sigma_8_sqr()    
{
  static double N_SPEC = -42.;
  static double res = -42.;
  static double OMEGA_M = -42.;
  static double H0 = -42.;
  static double OMB = -42.;
  double integral;
  
  if (OMEGA_M != cosmology.Omega_m ||H0 != cosmology.h0 || OMB!= cosmology.omb|| N_SPEC != cosmology.n_spec)
  {
    // full integral - method of choice!
    integral = int_GSL_integrate_qag(int_for_sigma_8,NULL,1e-4,1e6,NULL,2048);
    res = 4.5/constants.pi/constants.pi*integral;   //see PD97, eq. 29
    
    H0=cosmology.h0;
    OMB=cosmology.omb;
    OMEGA_M=cosmology.Omega_m;
    N_SPEC = cosmology.n_spec;
  }
  assert(res>0.0);
  return res;
}



/* ============================================================ *
 * Calculates k_NL, n_eff, n_cur.				*
 * ============================================================ */ 

double int_for_wint2_knl(double logk)
{
  double krsqr, k;
  
  k = exp(logk);
  krsqr = DSQR(k*global.rglob);
  //printf("Hallo\n");
  return Delta_L_tab(k)*exp(-krsqr);
}

double int_for_wint2_neff(double logk)
{
  double krsqr, k;
  
  k = exp(logk);
  krsqr = DSQR(k*global.rglob);
  return Delta_L_tab(k)*2.0*krsqr*exp(-krsqr);
}

double int_for_wint2_ncur(double logk)
{
  double krsqr, k;
  
  k = exp(logk);
  krsqr = DSQR(k*global.rglob);
  return Delta_L_tab(k)*4.0*krsqr*(1.0-krsqr)*exp(-krsqr);
}




void wint2(double r, double *sig, double *d1, double *d2, double amp, int onlysig)
{
  const double kmin = 1.e-2;
  double kmax, logkmin, logkmax, s1, s2, s3;
  
  /* choose upper integration limit to where filter function dropped
   * substantially */
  kmax  = sqrt(5.*log(10.))/r;
  if (kmax<8000.0) kmax = 8000.0;
  
  logkmin = log(kmin);
  logkmax = log(kmax);
  global.rglob = r;
  if (onlysig==1) {
    //printf("Hallo\n");
    s1   = sm2_qromb(int_for_wint2_knl, logkmin, logkmax);
    *sig = amp*sqrt(s1);
  } else s1 = DSQR(1/amp);   /* sigma = 1 */
    
    if (onlysig==0) {
      s2  = sm2_qromb(int_for_wint2_neff, logkmin, logkmax);
      s3  = sm2_qromb(int_for_wint2_ncur, logkmin, logkmax);
      *d1 = -s2/s1;
      *d2 = -DSQR(*d1) - s3/s1;
    }
}

void nonlinscale(double ampsqr,double *rknl,double *rneff,double *rncur,int *golinear)
{
  double d1,d2,sig,diff,logr1,logr2,rmid,logr1start,logr2start,logrmidtmp,logrmid,amp;  
  int iter; 
  const double logstep = 5.0;
  const int itermax  = 40;
  const double kNLstern = 1.e6;       /* h/Mpc */
  
  logr1 = -2.0;
  logr2 =  3.5;
  
  iterstart:
  
  logr1start = logr1;
  logr2start = logr2; 
  //printf("iterstart %le %le\n",logr1start,logr2start);
  
  iter=0;
  do {
    logrmid=(logr2+logr1)/2.0;
    rmid=pow(10,logrmid);
    //printf("test%le %le %le\n",rmid,logr2,logr1);
    amp=sqrt(ampsqr);
    //printf("Hallo1\n");
    wint2(rmid, &sig, 0x0, 0x0, amp, 1);
    diff=sig-1.0;	
    //printf("%d %le\n",iter,diff);  
    if (diff>0.001){
      logr1=log10(rmid);
      //printf("hallo%le\n",xlogr1);
    }
    if (diff<-0.001){ 
      logr2=log10(rmid);
      //printf("hallo2%le\n",xlogr2);
    }
  } while (fabs(diff)>=0.001 && ++iter<itermax);
  
  if (iter>=itermax) {
    logrmidtmp = (logr2start+logr1start)/2.0;
    //printf("%le %le\n",logrmid,logrmidtmp);
    if (logrmid<logrmidtmp) {
      logr1 = logr1start-logstep;
      logr2 = logrmid;
    } else if (logrmid>=logrmidtmp) {
      logr1 = logrmid;
      logr2 = logr2start+logstep;
    }
    /* non-linear scale far beyond maximum scale: set flag 
     *    golinear */
    if (1/pow(10, logr1)>kNLstern) {
      *golinear = 1;
      //goto after_wint;
    } else  {
      goto iterstart;
    }
  }
  if (*golinear ==0){
    wint2(rmid, &sig, &d1, &d2, amp, 0);
    //printf("Hallo\n");
    *rknl  = 1./rmid;
    *rneff = -3-d1;
    *rncur = -d2;
  }
  //after_wint:
  //printf("Hallo\n");
}                  

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//c halo model cosmology.nonlinear fitting formula as described in 
//c Appendix C of Smith et al. (2002)

void halofit(double rk,double rn,double rncur,double rknl,double plin,double om_m, double om_v,double *pnl,double aa)   //double aa neu!
{
  double gam,a,b,c,xmu,xnu,z0,beta,f1,f2,f3;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac,pq,ph;
  double nsqr;
  
  nsqr = rn*rn;
  gam = 0.86485 + 0.2989*rn + 0.1631*rncur;
  a = 1.4861 + 1.83693*rn + 1.67618*nsqr + 0.7940*rn*nsqr
  + 0.1670756*nsqr*nsqr - 0.620695*rncur;
  a = pow(10,a);
  b = pow(10,(0.9463+0.9466*rn+0.3084*nsqr-0.940*rncur));
  c = pow(10,(-0.2807+0.6669*rn+0.3214*nsqr-0.0793*rncur));
  xmu = pow(10,(-3.54419+0.19086*rn));
  xnu = pow(10,(0.95897+1.2857*rn));
  z0 = 1.38848+0.3701*rn-0.1452*nsqr;
  beta = 0.8291+0.9854*rn+0.3400*nsqr;
  
  if(fabs(1-om_m)>0.01)
  {
    //printf("non EDS om_m=%le\n",om_m);
    
    f1a = pow(om_m,(-0.0732));
    f2a = pow(om_m,(-0.1423));
    f3a = pow(om_m,(0.0725));
    f1b = pow(om_m,(-0.0307));
    f2b = pow(om_m,(-0.0585));
    f3b = pow(om_m,(0.0743));
    frac = om_v/(1.-om_m); 
    
    //new - interpolation in eq. of state parameters - treat as open universe
    double we,w_de;
    w_de=cosmology.w0+cosmology.wa*(1.-aa);   //w_z(z);
    we=(frac*w_de)+((1.0-frac)*(-1.0/3.0));
    frac=(-1.0)*((3.0*we)+1.0)/2.0;
    //end of new
    
    
    f1 = frac*f1b + (1-frac)*f1a;
    f2 = frac*f2b + (1-frac)*f2a;
    f3 = frac*f3b + (1-frac)*f3a;
  }
  else /* EdS Universe */
  {
    //printf("EDS om_m=%le\n",om_m);
    f1 = 1.0;
    f2 = 1.;
    f3 = 1.;
  }
  y = (rk/rknl);
  ph = a*pow(y,f1*3)/(1+b*pow(y,f2)+pow(f3*c*y,3-gam));
  ph = ph/(1+xmu/y+xnu/y/y);
  pq = plin*pow(1+plin,beta)/(1+plin*z0)*exp(-y/4.0-y*y/8.0);
  *pnl = pq + ph;
  
  assert(finite(*pnl));
}


/* slope in the highly cosmology.nonlinear regime, c.f. Smith et al (2002) eq. (61) 
 */
double slope_NL(double rn, double rncur, double om_m, double om_v)
{
  double gam, f1a, f1b, frac, f1;
  
  gam = 0.86485 + 0.2989*rn + 0.1631*rncur;
  if(fabs(1-om_m)>0.01) {
    f1a = pow(om_m,(-0.0732));
    f1b = pow(om_m,(-0.0307));
    frac = om_v/(1.-om_m);  
    f1 = frac*f1b + (1-frac)*f1a;
  } else {
    f1 = 1.0;
  }
  return 3.0*(f1-1.0) + gam - 3.0;
}      


void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da){     
  
  double plin,pnl,rk,rknl,rneff,rncur,om_m,om_v,amp,ampsqr,grow0,aa,klog,Pdelta;
  
  int i,j,golinear;
  grow0=growfac(1.,1.,1.);
  aa = limits.a_min;
  //  printf("New P_delta\n");
  for (i=0; i<table_N_a; i++, aa +=da) {
    if(aa>1.0) aa=1.0;
    omega_a(aa,&om_m,&om_v);
    amp=growfac(aa,1.,1.)/grow0;
    ampsqr=amp*amp;
    golinear=0;
    nonlinscale(ampsqr,&rknl,&rneff,&rncur,&golinear);
    //printf("New nonlin\n");
    
    klog = logkmin;
    for (j=0; j<table_N_k; j++, klog += dk) {
      rk=exp(klog);
      plin=ampsqr*Delta_L_tab(rk);
      if(golinear==0) {
	halofit(rk,rneff,rncur,rknl,plin,om_m,om_v,&pnl,aa);
	Pdelta=pnl;
      }
      else{
	Pdelta=plin;
      }
      //table_P_NL[i][j]=log(Pdelta*(1.0+2.0*(rk/10)*(rk/10))/(1.0+(rk/10)*(rk/10)));
      table_P_NL[i][j]=log(Pdelta);
    }
  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double Delta_nl_halo(double a,double k_NL)
{     
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double N_SPEC  = -42.;
  static double OMB   = -42.;
  static double H0  = -42.;
  static double SIGMA_8 = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double **table_P_NL = 0;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_P_NL_calibrate=0;
  
  double aa,klog,val;
  static double COSMO_min[5] = 
  {1.2000000000000000e-01,   2.1499999999999998e-02,   8.4999999999999998e-01,   6.1599999999999999e-01,  -1.3000000000000000e+00};
  static double COSMO_max[5] = {1.550000e-01,2.350000e-02,1.050000e+00,9.000000e-01,-7.000000e-01};
  
  double COSMO[7],COSMO_calib[6],stuff[4],kstar[1995],p_coyote[1995],scale_large_k,scale_small_k;
  int type=1;
  int i,j;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  ||W0 != cosmology.w0 || WA != cosmology.wa || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || SIGMA_8 != cosmology.sigma_8) {
    if (table_P_NL!=0) sm2_free_matrix(table_P_NL,0, table_N_a-1, 0, table_N_k-1);     
    if (table_P_NL_calibrate!=0) sm2_free_matrix(table_P_NL_calibrate,0, table_N_a-1, 0, table_N_k-1);     
    table_P_NL = sm2_matrix(0, table_N_a-1, 0, table_N_k-1);     
    table_P_NL_calibrate = sm2_matrix(0, table_N_a-1, 0, table_N_k-1);    
    
    da = (1. - limits.a_min)/(table_N_a-1.);
    logkmin = log(limits.Pdelta_halo_k_min);
    logkmax = log(limits.Pdelta_halo_k_max);
    dk = (logkmax - logkmin)/(table_N_k-1.);
    
    printf("%le %le %le %le %le %le %le\n",cosmology.Omega_m,cosmology.omb,cosmology.n_spec,cosmology.sigma_8,cosmology.w0,cosmology.wa,cosmology.h0);
    Delta_halofit(table_P_NL,logkmin, logkmax, dk, da);
    
    if (configuration.COYOTE_UNIVERSE_CALIBRATION==1){
      
      COSMO[0] = cosmology.Omega_m;
      COSMO[1] = cosmology.omb;
      COSMO[2] = cosmology.n_spec;
      COSMO[3] = cosmology.sigma_8;
      COSMO[4] = cosmology.w0;
      COSMO[5] = cosmology.wa;
      COSMO[6] = cosmology.h0;
      
      COSMO_calib[0]=0.1344;
      COSMO_calib[1]=0.02246;
      COSMO_calib[2]=cosmology.n_spec;
      COSMO_calib[3]=cosmology.sigma_8;
      COSMO_calib[4]= -1.000;
      
      if (COSMO[2] <= COSMO_min[2]) COSMO_calib[2]=COSMO_min[2]; 
      if (COSMO[3] <= COSMO_min[3]) COSMO_calib[3]=COSMO_min[3]; 
      if (COSMO[2] >= COSMO_max[2]) COSMO_calib[2]=COSMO_max[2];
      if (COSMO[3] >= COSMO_max[3]) COSMO_calib[3]=COSMO_max[3];
      cosmology.omh2=COSMO_calib[0];
      cosmology.ombh2=COSMO_calib[1];
      cosmology.n_spec=COSMO_calib[2];
      cosmology.sigma_8=COSMO_calib[3];	
      cosmology.w0=COSMO_calib[4];
      cosmology.wa=0.0;
      
      getH0fromCMB(COSMO_calib, stuff);
      if((stuff[3]>0.715) && (stuff[3]<0.703)) printf("error in coyote calibration scheme h0=%le\n",stuff[3]);
      cosmology.h0=stuff[3];
      cosmology.Omega_m=COSMO_calib[0]/cosmology.h0/cosmology.h0;
      cosmology.Omega_v=1-cosmology.Omega_m;
      cosmology.omb=COSMO_calib[1]/cosmology.h0/cosmology.h0;
      //printf("%le %le %le %le %le %le\n",COSMO_calib[0],COSMO_calib[1],COSMO_calib[2],COSMO_calib[3],COSMO_calib[4],COSMO_calib[5]);
      
      Delta_halofit(table_P_NL_calibrate,logkmin, logkmax, dk, da);
      aa = limits.a_min;
      for (i=0; i<table_N_a; i++, aa +=da) {
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *timspline = gsl_spline_alloc (gsl_interp_cspline, 1995);
	COSMO_calib[5]=(1.0/aa)-1.0;
	if(fabs(COSMO_calib[5])<1.e-10) {
	  COSMO_calib[5]=0.0;
	}
	if((aa >= a_min_coyote) && (aa <= a_max_coyote)){
	  emu(COSMO_calib,kstar,p_coyote, &type);
	  gsl_spline_init (timspline, kstar, p_coyote, 1995);
	  
	  klog=log(kstar[0]/cosmology.h0);
	  val=sm2_interpol2d(table_P_NL_calibrate, table_N_a, limits.a_min, 1., da, aa, table_N_k, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
	  scale_small_k=p_coyote[0]/exp(val);
	  
	  klog=log(kstar[1994]/cosmology.h0); 
	  val=sm2_interpol2d(table_P_NL_calibrate, table_N_a, limits.a_min, 1., da, aa, table_N_k, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
	  scale_large_k=p_coyote[1994]/exp(val);
	  
	  klog = logkmin; // log k in h/MPC
	  for (j=0; j<table_N_k; j++, klog += dk) {
	    if ((klog >= log(k_min_coyote/cosmology.h0)) && (klog <= log(k_max_coyote/cosmology.h0))){
	      table_P_NL[i][j]=log(gsl_spline_eval(timspline, exp(klog)*cosmology.h0, acc))-table_P_NL_calibrate[i][j]+table_P_NL[i][j];
	      
	      //		fprintf(F,"%le %le\n",exp(klog)*cosmology.h0,exp(table_P_NL[i][j]));
	    }
	    if(klog>log(k_max_coyote/cosmology.h0)){
	      table_P_NL[i][j]=table_P_NL[i][j]+log(scale_large_k);
	    }
	    if(klog<log(k_min_coyote/cosmology.h0)){
	      table_P_NL[i][j]=table_P_NL[i][j]+log(scale_small_k);
	    }
	  }
	}
	gsl_spline_free (timspline);
	gsl_interp_accel_free (acc);      
      }
      // reset cosmology  
      cosmology.Omega_m=COSMO[0];
      cosmology.Omega_v=1.0-cosmology.Omega_m;
      cosmology.omb=COSMO[1];
      cosmology.n_spec=COSMO[2];
      cosmology.sigma_8=COSMO[3];	
      cosmology.w0=COSMO[4];
      cosmology.wa=COSMO[5];
      cosmology.h0=COSMO[6];
      
    }
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    N_SPEC  = cosmology.n_spec;
    OMB = cosmology.omb;
    H0=cosmology.h0;
    SIGMA_8 = cosmology.sigma_8; 
    W0      = cosmology.w0;
    WA      = cosmology.wa;
    //		printf("P_delta_table finished\n");
  }    
  klog = log(k_NL); 
  //   printf("%le %le %le %d %d %le %le %le\n",a,limits.a_min,da,table_N_a,table_N_k, logkmin, logkmax, dk);
  val = sm2_interpol2d(table_P_NL, table_N_a, limits.a_min, 1., da, a, table_N_k, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
  return exp(val); 
  // returns the dimensionless power spectrum as a function of scale factor a and k in units of h/Mpc 
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*============================================================
 *see BS 2.41 bzw Logbook for detailed calculation of w from a.*/	

double int_for_chi(double a){
  double res,asqr;
  asqr=a*a;	
  res= 1./sqrt(a*cosmology.Omega_m + asqr*(1.-cosmology.Omega_m -cosmology.Omega_v ) + asqr*asqr*omv_vareos(a)); 
  return res;
}


/*for the calculation of chi we have to integrate from a(z2)=a up to a(z1)=1, which means todays expansion factor*/	
double chi(double a){
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W_0= -42.;
  static double W_A = -42.;
  static double *table;
  static double da = 0.;
  double aa;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v || W_0 != cosmology.w0   || W_A!= cosmology.wa  ){
    da = (1.-limits.a_min)/(table_N_a-1.);
    aa = limits.a_min;
    if (table!=0) sm2_free_vector(table, 0, table_N_a-1);
    table   = sm2_vector(0, table_N_a-1);
    for (i=0; i<table_N_a-1; i++, aa+=da) {
      table[i] = sm2_qromb(int_for_chi, aa, 1.);
    }
    table[table_N_a-1] =.0;
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    W_0 = cosmology.w0 ;
    W_A = cosmology.wa ;
  }
  return sm2_interpol(table, table_N_a, limits.a_min, 1., da, a, 1.0, 1.0); // comoving distance in c/H_0
}


/*===============================calculating the angular diameter distance f_K BS01 2.4, 2.30: f_K is a radial function that, depending on the curvature of the Universe, is a trigonometric, linear, or hyperbolic function of chi  */
double f_K(double chi)
{  
  double K, K_h, f;
  K = (cosmology.Omega_m   + cosmology.Omega_v  - 1.);
  if (K > precision.eps1) {           /* open */
    K_h = sqrt(K); // K in units H0/c see BS eq. 2.30
    f = 1./K_h*sin(K_h*chi);
    //printf("open\n");
  } else if (K < -precision.eps1) {   /* closed */
    K_h = sqrt(-K); 
    f = 1./K_h*sinh(K_h*chi);
    //printf("closed K=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  } else {                     /* flat */
    f = chi;
    //printf("flatK=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  }
  return f;
}

/*************************************************************************************/
/*      Redshift distribution with cuts, etc. (BG)                                     
 */
/************************************************************************************/
double p_zspec(double zs, double zp){
	static double init = -42.;
	static double **table_z = 0;
	static double zmin = 0.04, zmax = 1.96, dz = 0.08;
	double norm;
	if (init != 1.0){
		if (!table_z) {table_z = sm2_matrix(0, table_N_zs-1, 0, table_N_zp-1);}
		FILE *ein;
		int i,j;
		double a2,a3,a4,z,z1;
		ein = fopen("./distr_zphot_zspec","r");
		for (i = 0; i< table_N_zs; i++){
			for (j = 0; j< table_N_zp; j++){
				fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
				table_z[i][j] = a4;
			}
		}
		init = 1.0;
	}
	if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0;}
	else{ return sm2_interpol2d(table_z,table_N_zs, zmin, zmax,dz,zs,table_N_zp, zmin, zmax,dz,zp,1.0,1.0)/dz;}
}
double int_for_zdistr_photoz_inner(double z, void *params)
{
	double *array = (double*)params;
	double zz=z/redshiftshear.z0;
	return pow(zz,redshiftshear.alpha)*exp(-pow(zz,redshiftshear.beta_p))*p_zspec(array[2],z);
}

double int_for_zdistr_photoz(double ztrue, void *params)
{
	double *array = (double*)params;
    double result,error;
	array[2] = ztrue;
    gsl_integration_workspace *w;
    gsl_function H;
	
    w = gsl_integration_workspace_alloc (1000);
	
    H.function = &int_for_zdistr_photoz_inner;
    H.params = (void*)array;
    
    gsl_integration_qag (&H,array[0],array[1], 0, 1.e-3, 1000, GSL_INTEG_GAUSS41,
						 w, &result, &error); 
    gsl_integration_workspace_free(w);
    return result;
}
double int_for_zdistr_mock_photoz_inner(double z, void *params)
{
	int i; 
	double val,*z_vector,*Nz_vector,space1,space2,ztrue;
	double *array = (double*)params;
	char filename[200];
	FILE *ein;
	int static table=0;
	double static zhisto_max,zhisto_min;
	ztrue = array[2];
	if (table!=1){
		z_vector=sm2_vector(0, redshiftshear.histogram_zbins-1);
		Nz_vector=sm2_vector(0, redshiftshear.histogram_zbins-1);
		const gsl_interp_type *redshift_t5 = gsl_interp_cspline;
		zaccel5 = gsl_interp_accel_alloc();
		redshift_distrib_spline = gsl_spline_alloc (redshift_t5,redshiftshear.histogram_zbins);
		sprintf(filename,"%s",redshiftshear.REDSHIFT_FILE);	  
		printf("%s\n",filename);
		ein=fopen(filename,"r");
		
		for (i=0;i<redshiftshear.histogram_zbins;i++){
			fscanf(ein,"%le %le %le %le\n",&space1,&z_vector[i],&space2,&Nz_vector[i]);
			printf("%le %le\n",z_vector[i],Nz_vector[i]);
		}
		fclose(ein);
		table=1;
		gsl_spline_init (redshift_distrib_spline, z_vector,Nz_vector,redshiftshear.histogram_zbins);    
		zhisto_max=z_vector[redshiftshear.histogram_zbins-1];
		zhisto_min=z_vector[0];
		
		sm2_free_vector(z_vector,0,redshiftshear.histogram_zbins-1);    
		sm2_free_vector(Nz_vector,0,redshiftshear.histogram_zbins-1);
		
		printf("USING 4 column z-distrib\n");
	}
	if((z<zhisto_min) && (z>=redshiftshear.zdistrpar_zmin)) z=0.11;
	if((z>zhisto_max) && (z<=redshiftshear.zdistrpar_zmax)) z=1.99;   
	val= gsl_spline_eval(redshift_distrib_spline,z,zaccel5)*p_zspec(ztrue,z);
	//printf("%le\n",val);
	if (isnan(val))  printf("Redshift range exceeded: z=%le\n",z);
	return val;
}
double int_for_zdistr_mock_photoz(double ztrue, void *params)
{
	double *array = (double*)params;
    double result,error;
	array[2] = ztrue;
    gsl_integration_workspace *w;
    gsl_function H;

    w = gsl_integration_workspace_alloc (1000);
	
    H.function = &int_for_zdistr_mock_photoz_inner;
    H.params = (void*)array;
    
    gsl_integration_qag (&H,array[0],array[1], 0, 1.e-3, 1000, GSL_INTEG_GAUSS41,
						 w, &result, &error); 
    gsl_integration_workspace_free(w);
    return result;
}
double zdistr_photoz(double z,int i) //returns p(ztrue | i), works with binned or analytic distributions; i =-1 -> no tomography; i>= 0 -> tomography bin i
{
	static double norm = 0.;
	static double BETA_P = -42.;
	static double ALPHA  = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	double x, f,array[3];
	switch (i) {
		case 0:array[0] = redshiftshear.zdistrpar_zmin1;array[1] = redshiftshear.zdistrpar_zmax1;
			break;
		case 1:array[0] = redshiftshear.zdistrpar_zmin2;array[1] = redshiftshear.zdistrpar_zmax2;
			break;
		case 2:array[0] = redshiftshear.zdistrpar_zmin3;array[1] = redshiftshear.zdistrpar_zmax3;
			break;
		case 3:array[0] = redshiftshear.zdistrpar_zmin4;array[1] = redshiftshear.zdistrpar_zmax4;
			break;
		case 4:array[0] = redshiftshear.zdistrpar_zmin5;array[1] = redshiftshear.zdistrpar_zmax5;
			break;
		default:array[0] = redshiftshear.zdistrpar_zmin; array[1] = redshiftshear.zdistrpar_zmax;
			break;
	}
	//First, compute the normalization
	if (ALPHA != redshiftshear.alpha || BETA_P != redshiftshear.beta_p || Z0 != redshiftshear.z0 || ZMIN !=redshiftshear.zdistrpar_zmin || ZMAX !=redshiftshear.zdistrpar_zmax)
	{
		gsl_integration_workspace *w;
		gsl_function H;
		
		w = gsl_integration_workspace_alloc (1000);
		
		H.params = (void*)array;
		
		if(redshiftshear.histogram_zbins != 0 )
		{ 			
			H.function = &int_for_zdistr_mock_photoz;
		}
		
		if((redshiftshear.beta_p>0.) && (redshiftshear.histogram_zbins == 0 ))
		{
				H.function = &int_for_zdistr_photoz;
		}
		gsl_integration_qag (&H,redshiftshear.zdistrpar_zmin,redshiftshear.zdistrpar_zmax, 0, 1.e-3, 1000, GSL_INTEG_GAUSS41,w, &norm, &error); 
		gsl_integration_workspace_free(w);
		ALPHA  = redshiftshear.alpha;
		BETA_P = redshiftshear.beta_p;
		Z0     = redshiftshear.z0;
		ZMIN   = redshiftshear.zdistrpar_zmin;
		ZMAX   = redshiftshear.zdistrpar_zmax;
	}
	
	
	if(redshiftshear.histogram_zbins != 0)
	{
		if((redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax) && (z>redshiftshear.zdistrpar_zmax || z<redshiftshear.zdistrpar_zmin)) return 0.0;
		return int_for_zdistr_mock_photoz(z)/norm;
	}
	
	if((redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax) && (z>redshiftshear.zdistrpar_zmax || z<redshiftshear.zdistrpar_zmin)) return 0.0;
	return norm*int_for_zdistr_photoz(z);
}

double int_for_zdistr(double z)
{
	double zz=z/redshiftshear.z0;
	return pow(zz,redshiftshear.alpha)*exp(-pow(zz,redshiftshear.beta_p));
}


double int_for_zdistr_mock(double z)
{
  int i; 
  double val,*z_vector,*Nz_vector,space1,space2;
  char filename[200];
  FILE *ein;
  int static table=0;
  double static zmax,zmin;
  
  if (table!=1){
    z_vector=sm2_vector(0, redshiftshear.histogram_zbins-1);
    Nz_vector=sm2_vector(0, redshiftshear.histogram_zbins-1);
    const gsl_interp_type *redshift_t5 = gsl_interp_cspline;
    zaccel5 = gsl_interp_accel_alloc();
    redshift_distrib_spline = gsl_spline_alloc (redshift_t5,redshiftshear.histogram_zbins);
    sprintf(filename,"%s",redshiftshear.REDSHIFT_FILE);	  
    printf("%s\n",filename);
    ein=fopen(filename,"r");
    
    for (i=0;i<redshiftshear.histogram_zbins;i++){
      fscanf(ein,"%le %le %le %le\n",&space1,&z_vector[i],&space2,&Nz_vector[i]);
      printf("%le %le\n",z_vector[i],Nz_vector[i]);
    }
    fclose(ein);
    table=1;
    gsl_spline_init (redshift_distrib_spline, z_vector,Nz_vector,redshiftshear.histogram_zbins);    
    zmax=z_vector[redshiftshear.histogram_zbins-1];
    zmin=z_vector[0];
    
    sm2_free_vector(z_vector,0,redshiftshear.histogram_zbins-1);    
    sm2_free_vector(Nz_vector,0,redshiftshear.histogram_zbins-1);
    
    printf("USING 4 column z-distrib\n");
  }
  if((z<zmin) && (z>=redshiftshear.zdistrpar_zmin)) z=0.11;
  if((z>zmax) && (z<=redshiftshear.zdistrpar_zmax)) z=1.99;   
  val= gsl_spline_eval(redshift_distrib_spline,z,zaccel5);
  //printf("%le\n",val);
  if (isnan(val))  printf("Redshift range exceeded: z=%le\n",z);
  return val;
}


double zdistr(double z)
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

/*************************************************************************************/
/* ! global variable 'global.aglob' is needed here ! */
double int_for_g(double aprime)
{
  double chi_glob, chi_prime,val;
  chi_glob = chi(global.aglob);
  chi_prime = chi(aprime);
  val=zdistr(1./aprime-1.)*f_K(chi_prime-chi_glob)/f_K(chi_prime)/(aprime*aprime);
  
  return val;
}

/*redshift weighted lens efficiency factor for a source redshift distribution p_chi(d_chi) WL02(94)*/	
double g_source(double a)
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
  static double *table;
  
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
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W_0 != cosmology.w0   || W_A!= cosmology.wa  || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin) 
  {
    if (table!=0) sm2_free_vector(table, 0, table_N_a-1);
    table   = sm2_vector(0, table_N_a-1);
    da = (1.-limits.a_min)/(table_N_a-1.);
    table[0] = 0.0;
    aa = limits.a_min+da;
    /*case of redshift distribution */
    for (i=1;i<table_N_a-1;i++,aa+=da) {
      global.aglob = aa;
      table[i] = sm2_qromb(int_for_g, limits.a_min, aa);
    }
    table[table_N_a-1] = 1.;
    OMEGA_M = cosmology.Omega_m ;
    OMEGA_V = cosmology.Omega_v ;
    W_0 = cosmology.w0 ;
    W_A = cosmology.wa ;
    
    BETA_P  =   redshiftshear.beta_p;
    ALPHA   =   redshiftshear.alpha;
    Z0      =   redshiftshear.z0;
    Z_MAX   =   redshiftshear.zdistrpar_zmax;
    Z_MIN   =   redshiftshear.zdistrpar_zmin;
  }
  return sm2_interpol(table, table_N_a, limits.a_min, 1., da, a, 1.0, 1.0);
}


//============================================================	

double int_for_p_shear_shear(double a)
{
  double res,hoverh0,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  res= DSQR(g_source(a)/(a*a))/hoverh0;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; //k in units H0/c
  return res;
}


double P_shear_shear(double s)
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
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_vector(table, 0, table_N_s-1);
    table   = sm2_vector(0, table_N_s-1);
    for (i=0; i<table_N_s; i++, slog+=ds) {
      global.sglob = exp(slog);
      f1 = sm2_qromb(int_for_p_shear_shear, limits.a_min, 0.7); 
      f2 = sm2_qromo(int_for_p_shear_shear, .7, 1.0, sm2_midpnt);
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
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


/* ============================================================ *
 * Shear-shear correlation function xi_+ (pm=+1) and xi_- (pm=-1).   *
 * ============================================================ */


double xi_shear_shear(int pm, double theta)
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
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin || NONLINEAR != cosmology.nonlinear) 
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
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_shear(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}


void xi_via_hankel_shear_shear(double **xi, double *logthetamin, double *logthetamax)
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
    lP[i] = l*P_shear_shear(l);
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

double xi_shear_magnification(int pm, double theta)
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
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || SIGMA_8 != cosmology.sigma_8 || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 !=cosmology.h0 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax ||Z_MIN !=redshiftshear.zdistrpar_zmin || NONLINEAR != cosmology.nonlinear) 
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
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_shear_magnification(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return 2.0*res; //as P_shear_mag = 2 P_shear_shear
}


void xi_via_hankel_shear_magnification(double **xi, double *logthetamin, double *logthetamax)
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
    lP[i] = l*P_shear_shear(l);
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

double pf(double z)
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
      return 1.0; 
      
    }
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


double int_for_p_shear_position(double a)
{
  double res,ell, fK, k;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  
  res= g_source(a)*pf(1./a-1.)/a/fK;
  
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k; 
  return res;
}



double P_shear_position(double s)  //see Eq. 157 in Schneider 2006 WL
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
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2,fK,a, k;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr || NONLINEAR != cosmology.nonlinear)
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
	f1 = sm2_qromb(int_for_p_shear_position, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_shear_position, .7, 1.0, sm2_midpnt);
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
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
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


double xi_shear_position(int pm, double theta)
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


void xi_via_hankel_shear_position(double **xi, double *logthetamin, double *logthetamax)
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

double xi_magnification_position(int pm, double theta)
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


void xi_via_hankel_magnification_position(double **xi, double *logthetamin, double *logthetamax)
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
double int_for_p_position_position(double a)
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



double P_position_position(double s)  //see Eq. 157 in Schneider 2006 WL
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
	f1 = sm2_qromb(int_for_p_position_position, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_position_position, .7, 1.0, sm2_midpnt);
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



double xi_position_position(int pm, double theta)
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
    xi_via_hankel_position_position(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}



void xi_via_hankel_position_position(double **xi, double *logthetamin, double *logthetamax)
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
    lP[i] = l*P_position_position(l);
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


/* ============================================================ *
 * Convolution kernel for Hankel transform with bias q,		*
 * Bessel_mu. arg[0]=q, arg[1]=mu.				*
 * ============================================================ */

void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc)
{
  fftw_complex a1, a2, g1, g2;
  int           mu;
  double        mod, xln2, si, co, d1, d2, pref, q;
  q = arg[0];
  mu = (int)(arg[1]+0.1);
  
  /* arguments for complex gamma */
  a1[0] = 0.5*(1.0+mu+q);
  a2[0] = 0.5*(1.0+mu-q);
  a1[1] = 0.5*x; a2[1]=-a1[1];
  cdgamma(a1,&g1);
  cdgamma(a2,&g2);
  xln2 = x*constants.ln2;
  si   = sin(xln2);
  co   = cos(xln2);
  d1   = g1[0]*g2[0]+g1[1]*g2[1]; /* Re */
  d2   = g1[1]*g2[0]-g1[0]*g2[1]; /* Im */
  mod  = g2[0]*g2[0]+g2[1]*g2[1];
  pref = exp(constants.ln2*q)/mod;
  
  (*res)[0] = pref*(co*d1-si*d2);
  (*res)[1] = pref*(si*d1+co*d2);
}

/*=========================================================*/
/*============= IA routines ==============*/
/*=========================================================*/
double int_for_p_GI(double a)
{
  double res,ell,fK, k,D_a,GIz;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  //printf("hallo\n");
  D_a=growfac(a,1.,1.)/growfac(1.,1.,1.);
  //printf("%le %le\n",nuisance.oneplusz0_ia,nuisance.eta_ia);
  
  res= g_source(a)*pf(1./a-1.)/a/fK/D_a;
  //printf("%le %le\n",nuisance.oneplusz0_ia,nuisance.eta_ia);
  GIz=pow((1+(1./a-1))/(nuisance.oneplusz0_ia),nuisance.eta_ia);
  
  res=res*Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k*GIz; 
  return res;
}


double P_GI(double s)  //see Eq. 157 in Schneider 2006 WL
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
  static double NUISANCE =-42;
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr || NONLINEAR != cosmology.nonlinear|| NUISANCE!=nuisance
.eta_ia)
  {
    logsmin = log(limits.P_2_s_min);
    logsmax = log(limits.P_2_s_max);
    ds = (logsmax - logsmin)/(table_N_s - 1.);
    slog = logsmin;
    if (table!=0) sm2_free_vector(table, 0, table_N_s-1);
    table   = sm2_vector(0, table_N_s-1);
    
    for (i=0; i<table_N_s; i++, slog+=ds) {
      global.sglob = exp(slog);
 	f1 = sm2_qromb(int_for_p_GI, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_GI, .7, 1.0, sm2_midpnt);
	table[i]= log(cosmology.bias*cosmology.rcorr*2.0*constants.pi_sqr*(f1 + f2)); 
	//printf("calculating P_shear_position %d\n",i);
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
    Z_MAX   = redshiftshear.zdistrpar_zmax;
    Z_MIN   = redshiftshear.zdistrpar_zmin;
    BETA_P2  =   redshiftclustering.beta_p;
    ALPHA2   = redshiftclustering.alpha;
    Z02      =   redshiftclustering.z0;
    Z_MAX2   = redshiftclustering.zdistrpar_zmax;
    Z_MIN2   = redshiftclustering.zdistrpar_zmin;
    BIAS   = cosmology.bias;
    RCORR   = cosmology.rcorr;
    NUISANCE=nuisance.eta_ia;
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}


double xi_GI(int pm, double theta)
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
  static double NUISANCE =-42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 || BETA_P != redshiftshear.beta_p || ALPHA != redshiftshear.alpha || Z0 != redshiftshear.z0 ||Z_MAX !=redshiftshear.zdistrpar_zmax || Z_MIN !=redshiftshear.zdistrpar_zmin|| BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear || NUISANCE!=nuisance.eta_ia)
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
    NUISANCE=nuisance.eta_ia;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_GI(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[0], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  
  return res;
}


void xi_via_hankel_GI(double **xi, double *logthetamin, double *logthetamax)
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
    lP[i] = l*P_GI(l);
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


/*=========================================================*/
/*============= Intrinsic Alignment II-Term ==============*/
/*=========================================================*/


double int_for_p_II(double a)
{
  double res,hoverh0,ell, fK, k, IIz;
  if (a >= 1.0) sm2_error("a>=1 in int_for_p_2");
  
  ell       = global.sglob;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
  
  res= pf(1./a-1.)*pf(1./a-1.)*a*a*hoverh0/fK/fK;
  IIz=pow((1+(1./a-1))/(nuisance.oneplusz0_ia),2.0*nuisance.eta_ia);
  res= res*Delta_nl_halo(a,k/cosmology.coverH0)*IIz/k/k/k; 
  return res;
}




double P_II(double s) 
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
  static double NUISANCE=-42.;
  
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1,slog,f2,fK,a, hoverh0,k;
  int i;
  
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 ||  BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear || NUISANCE != nuisance.eta_ia)
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
	f1 = sm2_qromb(int_for_p_II, limits.a_min, 0.7); 
	f2 = sm2_qromo(int_for_p_II, .7, 1.0, sm2_midpnt);
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
    NUISANCE   = nuisance.eta_ia;
  }
  slog = log(s);
  f1 = sm2_interpol(table, table_N_s, logsmin, logsmax, ds, slog, cosmology.n_spec, cosmology.n_spec-4.0);
  return exp(f1);
}



double xi_II(int pm, double theta)
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
  static double NUISANCE   = -42.;
  
  static double NONLINEAR = -42;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res;
  if (OMEGA_M != cosmology.Omega_m   || OMEGA_V != cosmology.Omega_v  || W0 != cosmology.w0   || WA!= cosmology.wa  || N_SPEC != cosmology.n_spec || OMB !=cosmology.omb || H0 != cosmology.h0 || SIGMA_8 != cosmology.sigma_8 ||  BETA_P2 != redshiftclustering.beta_p || ALPHA2 != redshiftclustering.alpha || Z02 != redshiftclustering.z0 ||Z_MAX2 !=redshiftclustering.zdistrpar_zmax || Z_MIN2 !=redshiftclustering.zdistrpar_zmin|| BIAS != cosmology.bias ||RCORR !=cosmology.rcorr|| NONLINEAR != cosmology.nonlinear|| NUISANCE != nuisance.eta_ia)
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
    NUISANCE   = nuisance.eta_ia;
    
    if (table!=0) sm2_free_matrix(table, 0, 1, 0, table_N_thetaH-1);
    table   = sm2_matrix(0, 1, 0, table_N_thetaH-1);
    xi_via_hankel_II(table, &logthetamin, &logthetamax);
    dlogtheta = (logthetamax-logthetamin)/((double)table_N_thetaH-1.0);
  }
  res = sm2_interpol(table[(1-pm)/2], table_N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
  return res;
}



void xi_via_hankel_II(double **xi, double *logthetamin, double *logthetamax)
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
    lP[i] = l*P_II(l);
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


