#ifndef __COSMOLOGY_H
#define __COSMOLOGY_H

typedef struct {
  int transferfunction_EH99;
  int COYOTE_UNIVERSE_CALIBRATION;
  int TOMOGRAPHY;
  int IA;
  
}configpara;

typedef struct {
  double sglob;
  double aglob;
  double rglob;
  int reset_cov_tables;
}globalpara;

typedef struct {
  double A_ia;
  double beta_ia;
  double eta_ia;
  double loverl0_ia;
  double oneplusz0_ia;
  double c1rhocrit_ia;
}nuisancepara;

typedef struct {
  int shear_Nbin; // number of tomography bins
  int shear_Npowerspectra;// number of tomography power spectra=1+2+3+...+Nbin
  double shear_zmax[10]; 
  double shear_zmin[10];
  int clustering_Nbin; // number of tomography bins
  int clustering_Npowerspectra;// number of tomography power spectra=1+2+3+...+Nbin
  double clustering_zmax[10]; 
  double clustering_zmin[10];
}tomopara;

typedef struct {
  int shear_photoz;
  double shear_zdistrpar_zmin;
  double shear_zdistrpar_zmax;
  int shear_histogram_zbins;
  char shear_REDSHIFT_FILE[200];
  double shear_z0;		/* redshift distribution scale */
  double shear_alpha;  /* redshift distribution power index.*/
  double shear_beta_p;		/* redshift distribution power index. If beta_b = 0, a distribution with a single source redshift at z0 is assumed */
   
  int clustering_photoz;
  double clustering_zdistrpar_zmin;
  double clustering_zdistrpar_zmax;
  int clustering_histogram_zbins;
  char clustering_REDSHIFT_FILE[200];
  double clustering_z0;		
  double clustering_alpha;  /* redshift distribution power index.*/
  double clustering_beta_p;		/* redshift distribution power index. If beta_b = 0, a distribution with a single source redshift at z0 is assumed */
   
}redshiftpara;

typedef struct {
     double Omega_m;		/* matter density parameter                       */
     double Omega_v;		/* cosmogical constant parameter                  */
     double sigma_8;		/* power spectrum normalization                   */
     double n_spec;		/* spectral index of initial power spectrum       */
     int  nonlinear;         /* 0: linear power spectrum                       */
     double w0; //time dependent Dark energy parametrization zero order
     double wa; //time dependent Dark energy parametrization first order
     double omb; //Omega baryon
     double h0; //Hubble constant
     double coverH0; //units for comoving distances - speeds up code
     double omh2;
     double ombh2;
     double rho_m;      /* = 3 H_0^2/(8 pi G) Omega_m, mean comoving matter density */
     double bias; 
     double rcorr; 
     double f_NL; 
}cosmopara;

void set_cosmological_parameters_to_WMAP_7years_BAO_SN();
void set_cosmological_parameters_to_WMAP_7years_only();
void set_cosmological_parameters_to_DES_mocks();
void set_cosmological_parameters_to_tak();
void set_cosmological_parameters_to_SATO();
void set_cosmological_parameters_to_OWLS();
void set_cosmological_parameters_to_Great10();
void set_redshift_DES();
void set_redshift_DES_conti();
void set_redshift_CFHTLS();
void set_redshift_SATO();
void set_tomo_BCC();
void set_tomo_DES_conti();
void set_redshift_SDSS();
void set_redshift_BCC();

#endif
