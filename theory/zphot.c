#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "maths.h"
#include "zphot.h"
#include "cosmology.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

gsl_interp_accel *zaccel5;
gsl_spline *redshift_distrib_spline;

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
	5,
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
	0.0
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
	//  sprintf(redshiftshear.REDSHIFT_FILE,"BCC_z-distribution");
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




double p_zphot(double zp, double zs){
	static double init = -42.;
	static double **table_z = 0;
	static double zmin = 0.04, zmax = 1.96, dz = 0.08;
	double norm;
	if (init != 1.0){
		if (!table_z) {table_z = sm2_matrix(0, table_N_zp-1, 0, table_N_zs-1);}
		FILE *ein;
		int i,j;
		double a2,a3,a4,z,z1;
		ein = fopen("./distr_zphot_zspec","r");
		for (i = 0; i< table_N_zs; i++){
			for (j = 0; j< table_N_zp; j++){
				fscanf(ein,"%le %le %le %le %le",&z1,&a2,&a3,&a4,&z);
				table_z[j][i] = a3;
			}
		}
		init = 1.0;
	}
	if (zs < zmin || zp < zmin || zs> zmax ||zp> zmax){return 0;}
	else{	return sm2_interpol2d(table_z,table_N_zp, zmin, zmax,dz,zp,table_N_zs, zmin, zmax,dz,zs,1.0,1.0)/dz;}	
}
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
	double x, f,array[3],error;
	if (i < 0){array[0] = redshiftshear.zdistrpar_zmin; array[1] = redshiftshear.zdistrpar_zmax;}
	if (i>= 0 && i < sheartomo.Nbin){array[0] =sheartomo.zmin[i];array[1] = sheartomo.zmax[i];}
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
		return int_for_zdistr_mock_photoz(z,(void*)array)/norm;
	}
	
	if((redshiftshear.zdistrpar_zmin || redshiftshear.zdistrpar_zmax) && (z>redshiftshear.zdistrpar_zmax || z<redshiftshear.zdistrpar_zmin)) return 0.0;
	return int_for_zdistr_photoz(z,(void*)array)/norm;
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


int main(){
	double z;
	double a1,a2,a3,b1,b2,b3,t[3];
	int i;
	set_redshift_BCC();
	set_tomo_BCC();
	sprintf(redshiftshear.REDSHIFT_FILE,"../../z_histo_bin95_larger_z01");
	for (z =0.05; z<1.9;z+=.05){
		for (i = 0; i<=2; i++){
			if (z > sheartomo.zmin[i] && z< sheartomo.zmax[i]){t[i] = zdistr(z);}
			else{t[i]= 0.0;}
		}
		b1 = zdistr_photoz(z,0);
		b2 = zdistr_photoz(z,1);
		b3 = zdistr_photoz(z,2);
		printf("%e   %e %e  %e %e  %e %e\n",z,t[0],b1,t[1],b2,t[2],b3);
		a1 += b1*0.05;
		a2 += b2*0.05;
		a3 += b3*0.05;
	}
	printf("%e %e %e\n",a1,a2,a3);
}