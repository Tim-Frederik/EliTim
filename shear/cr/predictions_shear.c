#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "predictions_shear.h"
#include "../../theory/theory_all.h"
#include "../../theory/maths.h"
#include "../../theory/cosmology.h"
#include "../../theory/tomo.h"
#include "../../theory/survey.h"
#include "../../theory/halo.h"


extern con constants;
extern pre precision;
extern lim limits;
extern cosmopara cosmology;
extern sheartomopara sheartomo;
extern configpara configuration;
extern redshiftshearpara redshiftshear;
extern globalpara global;
extern nuisancepara nuisance;
extern sur survey;

#define ellmax 200000


void make_power_single_2D(char *OUT_FILE, double *ell, int Nbins)
{
  int j;
  char filename[300];
  FILE *aus;
  sprintf(filename,"%sp_shear_shear2D",OUT_FILE);
  aus=fopen(filename,"w");     
  for (j=0; j<Nbins; j++) {
    fprintf(aus,"%le %le\n",ell[j],P_shear_shear(ell[j]));
    //	printf("%le %le\n",ell[j],,P_shear_shear(ell[j]));
  }
  fclose(aus);
}

void make_power_single_tomo(char *OUT_FILE, double *ell, int Nbins)
{	
  int i,j,l;
  char filename[300];
  FILE *aus;
  for (i=0;i<sheartomo.Nbin; i++){
    redshiftshear.zdistrpar_zmin = sheartomo.zmin[i]; 
    redshiftshear.zdistrpar_zmax = sheartomo.zmax[i];					 
    for (j=i;j<sheartomo.Nbin; j++){
      redshiftshear.zdistrpar_zmin2 = sheartomo.zmin[j]; 
      redshiftshear.zdistrpar_zmax2 = sheartomo.zmax[j];			
      sprintf(filename,"%sp_shear_shear_tomo_z%d%d",OUT_FILE,i+1,j+1);
      aus=fopen(filename,"w");     
      for (l=0; l<Nbins; l++) {
	fprintf(aus,"%le %le\n",ell[l],P_shear_shear_tomo(ell[l]));
	//	printf("%le %le\n",ell[j],,P_shear_shear_tomo(ell[j]));
      }
      fclose(aus);
    }
  }    
}


void make_power_list_2D(char *OUT_FILE, double *ell,int Nbins, int START,int END, int N_para, char *PARA_FILE, int set_cosmo)
{
  
  int i,j;
  char filename[200];
  FILE *aus,*ein;
  
  double *Om,*Ob,*ns,*s8,*w0;
  
  Om= sm2_vector(0, N_para-1);
  Ob= sm2_vector(0, N_para-1);
  ns= sm2_vector(0, N_para-1);
  s8 = sm2_vector(0, N_para-1);
  w0= sm2_vector(0, N_para-1);
  
  ein=fopen(PARA_FILE,"r");
  for (i=0;i<N_para;i++){
    fscanf(ein,"%le %le %le %le %le\n",&Om[i],&Ob[i],&ns[i],&s8[i],&w0[i]);
  }
  fclose(ein);
  
  sprintf(filename,"%sp_shear_shear_2D_list%d_end%d",OUT_FILE,START,END);
  aus=fopen(filename,"w");
  
  for (j=START; j<END;j++){
    printf("Omega_m=%le sigma_8=%le\n",Om[j],s8[j]);
    set_cosmology(set_cosmo);
    cosmology.Omega_m=Om[j];
    cosmology.omb=Ob[j];
    cosmology.n_spec=ns[j];
    cosmology.sigma_8=s8[j];	
    cosmology.w0=w0[j];
    for (i=0; i<Nbins; i++) {
      fprintf(aus,"%le %le\n",ell[i],P_shear_shear(ell[i]));
      //	printf("%le %le\n",ell[j],,P_shear_shear(ell[i]));
    }
    printf("Power spectra Shear Shear Finished STEP %d of %d\n",i+1,END-START);
  }
  fclose(aus);
  
}

void make_power_list_tomo(char *OUT_FILE, double *ell,int Nbins, int START,int END, int N_para,char *PARA_FILE,int set_cosmo)
{
  
  int i,j,k,l;
  char filename[200];
  FILE *aus,*ein;
  
  double *Om,*Ob,*ns,*s8,*w0;
  
  Om= sm2_vector(0, N_para-1);
  Ob= sm2_vector(0, N_para-1);
  ns= sm2_vector(0, N_para-1);
  s8 = sm2_vector(0, N_para-1);
  w0= sm2_vector(0, N_para-1);
  
  ein=fopen(PARA_FILE,"r");
  for (i=0;i<N_para;i++){
    fscanf(ein,"%le %le %le %le %le\n",&Om[i],&Ob[i],&ns[i],&s8[i],&w0[i]);
  }
  fclose(ein);
  
  sprintf(filename,"%sp_shear_shear_tomo_list%d_end%d",OUT_FILE,START,END);
  aus=fopen(filename,"w");
  
  for (j=START; j<END;j++){
    printf("Omega_m=%le sigma_8=%le\n",Om[j],s8[j]);
    set_cosmology(set_cosmo);
    cosmology.Omega_m=Om[j];
    cosmology.omb=Ob[j];
    cosmology.n_spec=ns[j];
    cosmology.sigma_8=s8[j];	
    cosmology.w0=w0[j];
    for (i=0;i<sheartomo.Nbin; i++){
      redshiftshear.zdistrpar_zmin = sheartomo.zmin[i]; 
      redshiftshear.zdistrpar_zmax = sheartomo.zmax[i];					 
      for (k=i;k<sheartomo.Nbin; k++){
	redshiftshear.zdistrpar_zmin2 = sheartomo.zmin[k]; 
	redshiftshear.zdistrpar_zmax2 = sheartomo.zmax[k];			
	sprintf(filename,"%sp_shear_shear_tomo_z%d%d_list%d_%d",OUT_FILE,i+1,k+1,START,END);
	aus=fopen(filename,"w");     
	for (l=0; l<Nbins; l++) {
	  fprintf(aus,"%le %le\n",ell[l],P_shear_shear_tomo(ell[l]));
	  //	printf("%le %le\n",ell[j],,P_shear_shear_tomo(ell[j]));
	}
	fclose(aus);
      }
    }
  }
}  

void make_cov_power_single_2D(char *OUT_FILE,double *ell,double *dell,int Ncl)
{
  char filename[200];
  double prefac,A,n,noise_ss,a,b,c,noise,fsky,prefac3,prefac2,lin,quad,tl,th1,th2,thsv,Pi;
  int i,j;
  FILE *F1;
  
  A=survey.area*survey.area_conversion_factor;
  n=survey.n_gal*survey.n_gal_conversion_factor;
  noise_ss= survey.sigma_e*survey.sigma_e/2.0/n;
  
  fsky = survey.area/41253.0;
  
  prefac=2.0*constants.pi/A;
  prefac3= 3.0/2.0*cosmology.Omega_m;
  
  printf("surveyarea=%le sigma_e=%le n_gal=%le\n",survey.area,survey.sigma_e,survey.n_gal);
  printf("A=%le sigma_e_sqr/2/n=%le\n",A,noise_ss);
  
  sprintf(filename,"cov_new/cov_2D");
  F1=fopen(filename,"w");
  for(i=0; i<Ncl; i++){
    prefac2=prefac/ell[i]/dell[i];
    Pi=P_shear_shear(ell[i]);
    noise=prefac2*noise_ss*noise_ss;
    lin=prefac2*2.0*Pi*noise_ss;
    quad=prefac2*Pi*Pi;
    for(j=i; j<Ncl; j++){
	a=b=c=0.0;
	if(i==j){
	  a=noise;
	  b=lin;
	  c=quad;
	}
	printf("%le %le %le\n",a,b,c);
	tl = pow(prefac3,4.0)*project_tri_lin_cov(ell[i],ell[j],4)/(2.0*constants.twopi*fsky);
	printf("%le\n",tl);
	th1 = pow(prefac3,4.0)*project_tri_1h_cov(ell[i],ell[j],4)/(2.0*constants.twopi*fsky);
	printf("%le\n",th1);
	th2 = pow(prefac3,4.0)*project_tri_2h_cov(ell[i],ell[j],4)/(2.0*constants.twopi*fsky);
	printf("%le\n",th2);
	thsv=HSV_cov(ell[i],ell[j],fsky);
	fprintf(F1,"%le %le %le %le %le %le %le %le %le %le\n",ell[i],ell[j],a,b,c,tl,th1,th2,thsv,a+b+c+tl+th1+th2+thsv);
	printf("%le %le %le %le %le %le %le %le %le %le\n",ell[i],ell[j],a,b,c,tl,th1,th2,thsv,a+b+c+tl+th1+th2+thsv);
    }
  }
  fclose(F1);
}


void set_cosmology(int set_cosmo)
{
  if(set_cosmo==1){
    set_cosmological_parameters_to_WMAP_7years_BAO_SN();
  }
  if(set_cosmo==2){
    set_cosmological_parameters_to_WMAP_7years_only();
  }
  if(set_cosmo==3){
    set_cosmological_parameters_to_DES_mocks();
  }
  if(set_cosmo==4){
    set_cosmological_parameters_to_Great10();
  }
  if(set_cosmo==5){
    set_cosmological_parameters_to_OWLS();
  }
}


void set_redshift_distri(int set_redshift)
{
  if(set_redshift==1){
    set_redshift_CFHTLS();
  }
  if(set_redshift==2){
    set_redshift_SATO();
  }
  if(set_redshift==3){
    set_redshift_DES_conti(); //CFHTLS like but mean redshift =0.75 
  }
  if(set_redshift==4){
    set_redshift_SDSS(); //CFHTLS like but mean redshift =0.75 
  }
}

void write_parameters(int SINGLE, int COVARIANCES, char *PARA_FILE, char *OUT_FILE,int N_para,int START,int END,int Ncl, double lmax, double lmin, int set_redshift)
{
  printf("---------------------------------------\n");
  printf("---------RUN MODE PARAMETERS-----------\n");
  printf("\n");
  if (SINGLE==1)printf("Single cosmology mode\n");
  if (SINGLE==0)printf("List cosmology mode\n");
  if (COVARIANCES==0)printf("Prediction mode\n");
  if (COVARIANCES==1)printf("Covariance mode\n");
  
  printf("PARA_FILE: %s\n",PARA_FILE);
  printf("OUT_FILE: %s\n",OUT_FILE);
  printf("TOTAL # cosmologies: %d\n",N_para);
  printf("START: %d\n",START);
  printf("END: %d\n",END);
  printf("\n");
  printf("# bins CL: %d\n",Ncl);
  
  printf("\n");
  printf("lmax: %le\n",lmax);
  printf("lmin: %le\n",lmin);
  
  printf("---------------------------------------\n");
  printf("---------REDSHIFT PARAMETERS-----------\n");
  printf("\n");
  
  if(configuration.TOMOGRAPHY==0) {
    printf("Tomography: NO\n");
    printf("z_min: %le\n",redshiftshear.zdistrpar_zmin);
    printf("z_max: %le\n",redshiftshear.zdistrpar_zmax);
  }
  if(configuration.TOMOGRAPHY==1) {
    printf("Tomography: YES\n");
    printf("z1: %le - %le\n",sheartomo.zmin[0],sheartomo.zmax[0]);
    printf("z2: %le - %le\n",sheartomo.zmin[1],sheartomo.zmax[1]);
    printf("z3: %le - %le\n",sheartomo.zmin[2],sheartomo.zmax[2]);
    printf("z4: %le - %le\n",sheartomo.zmin[3],sheartomo.zmax[3]);
    printf("z5: %le - %le\n",sheartomo.zmin[4],sheartomo.zmax[4]);
  }
  
  if(redshiftshear.histogram_zbins==0){ 
    printf("Using Redshift Parametrization\n");
    printf("shear z0=%le\n",redshiftshear.z0);
    printf("shear beta: %le\n",redshiftshear.beta_p);
    printf("shear alpha: %le\n", redshiftshear.alpha);
  }
  if(redshiftshear.histogram_zbins!=0){ 
    printf("USING SHEAR REDSHIFT_FILE: %s\n",redshiftshear.REDSHIFT_FILE);
    printf("N_bins in dn/dz: %d\n",redshiftshear.histogram_zbins);
  }
}    


int main(int argc, char** argv)
{
  gsl_set_error_handler_off ();
  int i;
  char PARA_FILE[200],OUT_FILE[200];
  
  // RUN MODE setup
  int SINGLE=atoi(argv[1]);
  int COVARIANCES=atoi(argv[2]); 
  sprintf(PARA_FILE,"%s",argv[3]);
  sprintf(OUT_FILE,"%s",argv[4]);
  int N_para=atoi(argv[5]);
  int START=atoi(argv[6]);
  int END=atoi(argv[7]);
  
  int Ncl=atoi(argv[8]);
  double lmin= atof(argv[9]);
  double lmax= atof(argv[10]);
  
  // COSMOLOGY CODE setup
  configuration.transferfunction_EH99=atoi(argv[11]);
  configuration.COYOTE_UNIVERSE_CALIBRATION=atoi(argv[12]);
  int set_cosmo=atoi(argv[13]);
  set_cosmology(set_cosmo);
  
  // REDSHIFT setup
  int set_redshift=atoi(argv[14]);
  set_redshift_distri(set_redshift); //must be set to inittialize non-single plane mode
  redshiftshear.zdistrpar_zmin=atof(argv[15]);
  redshiftshear.zdistrpar_zmax=atof(argv[16]);
  redshiftshear.histogram_zbins=atoi(argv[17]);
  sprintf(redshiftshear.REDSHIFT_FILE,"%s",argv[18]);
    
  configuration.TOMOGRAPHY=atoi(argv[19]); 
  if(configuration.TOMOGRAPHY==1){
    sheartomo.Nbin=atoi(argv[20]);
    for (i=0;i<sheartomo.Nbin; i++){
      sheartomo.zmin[i]=atof(argv[21+i]);
      sheartomo.zmax[i]=atof(argv[22+i]);
    }
  }
  write_parameters(SINGLE,COVARIANCES,PARA_FILE,OUT_FILE,N_para,START,END,Ncl,lmax,lmin,set_redshift);
  
  //binning  calculation
  double logdl=(log(lmax)-log(lmin))/Ncl;
  double *ell, *dell;
  ell=sm2_vector(0,Ncl-1);
  dell=sm2_vector(0,Ncl-1);
  for(i=0; i<Ncl ; i++){
    ell[i]=exp(log(lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(lmin)+((i+1)*logdl)) - exp(log(lmin)+(i*logdl));
  }
  if((configuration.TOMOGRAPHY==0) && (SINGLE==1)) make_power_single_2D(OUT_FILE,ell,Ncl);     
  if((configuration.TOMOGRAPHY==0) && (SINGLE==0)) make_power_list_2D(OUT_FILE,ell,Ncl,START,END,N_para,PARA_FILE,set_cosmo);     
  if((configuration.TOMOGRAPHY==1) && (SINGLE==1)) make_power_single_tomo(OUT_FILE,ell,Ncl);     
  if((configuration.TOMOGRAPHY==1) && (SINGLE==0)) make_power_list_tomo(OUT_FILE,ell,Ncl,START,END,N_para,PARA_FILE,set_cosmo);
  
  set_survey_parameters_to_DES();
  if((COVARIANCES==1) && (configuration.TOMOGRAPHY==0) && (SINGLE==1)) make_cov_power_single_2D(OUT_FILE,ell,dell,Ncl);
  //   if((COVARIANCES==1) && (configuration.TOMOGRAPHY==0) && (SINGLE==0)) make_cov_list_2D(OUT_FILE,ell,dell,Ncl,,START,END,N_para,PARA_FILE,set_cosmo);
  //   if((COVARIANCES==1) && (configuration.TOMOGRAPHY==1) && (SINGLE==1)) make_cov_single_tomo(OUT_FILE,ell,dell, Ncl);
  //   if((COVARIANCES==1) && (configuration.TOMOGRAPHY==1) && (SINGLE==0)) make_cov_list_tomo(OUT_FILE,ell,dell,Ncl,START,END,N_para,PARA_FILE,set_cosmo); 
  
  return 0;
}


