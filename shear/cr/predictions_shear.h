#ifndef __SHEAR_H
#define __SHEAR_H

void make_power_single_2D(char *OUTFILE, double *ell, int Nbins);
void make_power_single_tomo(char *OUTFILE, double *ell, int Nbins);

void make_power_list_tomo(char *OUT_FILE, double *ell,int Nbins, int START,int END, int N_para,char *PARA_FILE,int set_cosmo);
void make_power_list_2D(char *OUT_FILE, double *ell,int Nbins, int START,int END, int N_para, char *PARA_FILE, int set_cosmo);
void make_cov_power_single_2D(char *OUT_FILE,double *ell,double *dell,int Ncl);
void set_cosmology(int set_cosmo);
void set_redshift_distri(int set_redshift);
void write_parameters(int SINGLE, int COVARIANCES, char *PARA_FILE, char *OUT_FILE,int N_para,int START,int END,int Ncl, double lmax, double lmin,int set_redshift);

#endif

