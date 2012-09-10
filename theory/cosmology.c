#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include "cosmology.h"

cosmopara cosmology = {
  0.0,
  0.0,
  0.0,
  0.0,
  0,
  0.0,
  0.0,
  0.0,
  0.0,
  2997.92458,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

redshiftshearpara redshiftshear = {
  1.0,
  0.0,
  0.0,
  0.0,
  3.0,
  "",
  0,
  0.0,
  3.0
};

sheartomopara sheartomo = {
  0,
  0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

redshiftclusteringpara redshiftclustering = {
  1.0,
  0.0,
  0.0,
  0.0,
  3.0,
  "",
  0
};

configpara configuration = {
  1,
  1,
  1
};


globalpara global = {
  0.0,
  0.0,
  0.0,
  -1
};

nuisancepara nuisance = {
  0.0, //A IA see Joachimi2012
  0.0, //beta IA see Joachimi2012
  0.0,  //eta_other IA see Joachimi2012
  0.96, //loverl0-ia
  1.54, //oneplusz0-ia MegaZ
  0.0134
};


void set_cosmological_parameters_to_WMAP_7years_BAO_SN()
{
  cosmology.Omega_m   = 0.272;
  cosmology.Omega_v   = 0.728;
  cosmology.sigma_8   = 0.807;
  cosmology.n_spec    = 0.961;
  
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0454463597;
  cosmology.h0=0.703;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  printf("Cosmology set to WMAP_7years_BAO_SN\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0; 
  cosmology.f_NL = 0.0;
}


void set_cosmological_parameters_to_WMAP_7years_only()
{
  cosmology.Omega_m     = 0.265;
  cosmology.Omega_v    = 0.735;
  cosmology.sigma_8   = 0.801;
  cosmology.n_spec    = 0.963;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.045;
  cosmology.h0=0.71;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  printf("Cosmology set to WMAP_7years_only\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0; 
  cosmology.f_NL = 0.0;
}

void set_cosmological_parameters_to_tak()
{
  cosmology.Omega_m     = 0.25;
  cosmology.Omega_v    = 0.75;
  cosmology.sigma_8   = 0.8;
  cosmology.n_spec    = 1.0;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.04;
  cosmology.h0=0.7;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  //  printf("Cosmology set to WMAP_7years_only\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0;
}



void set_cosmological_parameters_to_OWLS()
{
  cosmology.Omega_m   = 0.238;
  cosmology.Omega_v   = 0.762;
  cosmology.sigma_8   = 0.74;
  cosmology.n_spec    = 0.951;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0418;
  cosmology.h0=0.73;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  //printf("Cosmology set to OWLS\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0;
}

void set_cosmological_parameters_to_Great10()
{
  cosmology.Omega_m     = 0.25;
  cosmology.Omega_v    = 0.75;
  cosmology.sigma_8   = 0.8;
  cosmology.n_spec    = 0.963;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.0447927;
  cosmology.h0=0.8;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  printf("Cosmology set to GREAT 10\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0;
}

void set_cosmological_parameters_to_DES_mocks()
{
  cosmology.Omega_m   = 0.25;
  cosmology.Omega_v   = 0.75;
  cosmology.sigma_8   = 0.8;
  cosmology.n_spec    = 1.0;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.04;
  cosmology.h0=0.7;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;	
  printf("Cosmology set to DES_mocks\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0;
}


void set_cosmological_parameters_to_SATO()
{
  cosmology.Omega_m   = 0.238;
  cosmology.Omega_v   = 0.762;
  cosmology.sigma_8   = 0.76;
  cosmology.n_spec    = 0.958;
  cosmology.nonlinear = 1;
  cosmology.w0=-1.;
  cosmology.wa=0.;
  cosmology.omb=0.042;
  cosmology.h0=0.732;
  cosmology.omh2=cosmology.Omega_m*cosmology.h0*cosmology.h0;
  cosmology.ombh2=cosmology.omb*cosmology.h0*cosmology.h0;
  cosmology.rho_m = 7.4775e+21*cosmology.Omega_m;//cosmology.h0;
  printf("Cosmology set  to Sato\n");
  cosmology.bias=1.0;
  cosmology.rcorr = 1.0;
}


void set_redshift_DES()
{
  redshiftshear.z0        = 1.0;
  redshiftshear.beta_p    = 0.;
  redshiftshear.alpha   =    0.;
  redshiftshear.zdistrpar_zmin = 0.;
  redshiftshear.zdistrpar_zmax = 3.0;
  sprintf(redshiftshear.REDSHIFT_FILE,"dummy");
  redshiftshear.histogram_zbins=0;
  redshiftshear.zdistrpar_zmin2 = 0.;
  redshiftshear.zdistrpar_zmax2 = 3.0;
  
  redshiftclustering.z0        = 1.0;
  redshiftclustering.beta_p    = 0.0;
  redshiftclustering.alpha   =    0.0;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
  
}


void set_redshift_CFHTLS()
{
  redshiftshear.z0        = 0.555;
  redshiftshear.beta_p    = 1.197;
  redshiftshear.alpha   =   1.193;
  redshiftshear.zdistrpar_zmin = 0.;
  redshiftshear.zdistrpar_zmax = 3.0;
  sprintf(redshiftshear.REDSHIFT_FILE,"dummy");
  redshiftshear.histogram_zbins=0;
  redshiftshear.zdistrpar_zmin2 = 0.;
  redshiftshear.zdistrpar_zmax2 = 3.0;
  
  redshiftclustering.z0        = 0.555;
  redshiftclustering.beta_p    = 1.197;
  redshiftclustering.alpha   =   1.193;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
  
}


void set_redshift_DES_conti()
{
  redshiftshear.z0        = 0.47;
  redshiftshear.beta_p    = 1.197;
  redshiftshear.alpha   =   1.193;
  redshiftshear.zdistrpar_zmin = 0.;
  redshiftshear.zdistrpar_zmax = 3.0;
  sprintf(redshiftshear.REDSHIFT_FILE,"dummy");
  redshiftshear.histogram_zbins=0;
  redshiftshear.zdistrpar_zmin2 = 0.;
  redshiftshear.zdistrpar_zmax2 = 3.0;
  
  redshiftclustering.z0        = 0.47;
  redshiftclustering.beta_p    = 1.197;
  redshiftclustering.alpha   =   1.193;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
  printf("Redshift set to DES conti\n");
}

void set_redshift_BCC()
{
  redshiftshear.z0        = 0.47;
  redshiftshear.beta_p    = 1.197;
  redshiftshear.alpha   =   1.193;
  redshiftshear.zdistrpar_zmin = 0.1;
  redshiftshear.zdistrpar_zmax = 2.0;
  sprintf(redshiftshear.REDSHIFT_FILE,"dummy");
  redshiftshear.histogram_zbins=95;
  redshiftshear.zdistrpar_zmin2 = 0.1;
  redshiftshear.zdistrpar_zmax2 = 2.0;
  
  redshiftclustering.z0        = 0.47;
  redshiftclustering.beta_p    = 1.197;
  redshiftclustering.alpha   =   1.193;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
  printf("Redshift set to DES BCC\n");
}


void set_redshift_SATO()
{
  redshiftshear.z0        = 0.6;
  redshiftshear.beta_p    = 0.;
  redshiftshear.alpha   =    0.;
  redshiftshear.zdistrpar_zmin = 0.;
  redshiftshear.zdistrpar_zmax = 3.0;
  sprintf(redshiftshear.REDSHIFT_FILE,"dummy");
  redshiftshear.histogram_zbins=0;
}

void set_redshift_SDSS()
{
  redshiftshear.z0        = 0.6;
  redshiftshear.beta_p    = 0.0;
  redshiftshear.alpha   =   0.0;
  redshiftshear.zdistrpar_zmin = 0.0394737;
  redshiftshear.zdistrpar_zmax = 1.46053;
  sprintf(redshiftshear.REDSHIFT_FILE,"../nz-cunha_SDSS.dat");
  redshiftshear.histogram_zbins=19;
  redshiftshear.zdistrpar_zmin2 = 0.0394737;
  redshiftshear.zdistrpar_zmax2 = 1.46053;
  
  redshiftclustering.z0        = 0.6;
  redshiftclustering.beta_p    = 0.0;
  redshiftclustering.alpha   =   0.0;
  redshiftclustering.zdistrpar_zmin =0.0394737; 
  redshiftclustering.zdistrpar_zmax = 1.46053;
  sprintf(redshiftclustering.REDSHIFT_FILE,"../nz-cunha_SDSS.dat");
  redshiftclustering.histogram_zbins=19;
  printf("Redshift set to SDSS\n");
}


void set_tomo_BCC()
{
  sprintf(redshiftshear.REDSHIFT_FILE,"BCC_z-distribution");
  sheartomo.Nbin        = 3;
  sheartomo.zmax[0]      = 0.34;
  sheartomo.zmax[1]      = 0.9;
  sheartomo.zmax[2]      = 1.99;
  sheartomo.zmax[3]      = 0.0;
  sheartomo.zmax[4]      = 0.0;
  
  sheartomo.zmin[0]      = 0.01;
  sheartomo.zmin[1]      = 0.34;
  sheartomo.zmin[2]      = 0.9;
  sheartomo.zmin[3]      = 0.0;
  sheartomo.zmin[4]      = 0.0;
}

void set_tomo_DES_conti()
{
printf("set_tomo_DES_conti\n");
  sheartomo.Nbin        = 5;
  sheartomo.zmax[0]      = 1.06;
  sheartomo.zmax[1]      = 1.49;
  sheartomo.zmax[2]      = 1.83;
  sheartomo.zmax[3]      = 2.13;
  sheartomo.zmax[4]      = 3.0;
  
  sheartomo.zmin[0]      = 0.0;
  sheartomo.zmin[1]      = 1.06;
  sheartomo.zmin[2]      = 1.49;
  sheartomo.zmin[3]      = 1.83;
  sheartomo.zmin[4]      = 2.13;
}



















