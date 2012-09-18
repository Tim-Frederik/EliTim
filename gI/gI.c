#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "../theory/theory_all.h"
#include "../theory/maths.h"
#include "../theory/cosmology.h"
#include "../theory/halo.h"
#include "../theory/survey.h"
#include "../theory/tomo.h"

extern con constants;
extern pre precision;
extern lim limits;
extern cosmopara cosmology;
extern configpara configuration;
extern redshiftpara redshift;
extern tomopara tomo;
extern globalpara global;
extern nuisancepara nuisance;
extern sur survey;


double P_delta (double a, double k){
	
	return Delta_nl_halo(a,k/cosmology.coverH0)/k/k/k;
}
double int_for_C_gg_photoz (double a, void *params){
	double res,k,fK,D_a,GIz,hoverh0;
	double *array = (double*)params;
	
	fK     = f_K(chi(a));
	k      = array[0]/fK;
	D_a=growfac(a,1.,1.)/growfac(1.,1.,1.);
	hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
	
	res= a*a*hoverh0/fK/fK;
	
	res = res*P_delta(a, k);
	if (array[1] > 0.5){res = res*pow(p_zspec(1./a-1.,array[3])/(a*a),2.0);}
	return res;
}

double C_gg_photoz(double l, double zm, double DZ, double zphot){
	double RES,res,k,fK,D_a,GIz,hoverh0,z,a,norm,dz;
	double array[3] = {l,zphot,zm};
	dz = 0.02;
	RES = 0;
	norm = 0.;
	for (z = zm-DZ/2.0+dz/2.0; z< zm+DZ/2.0; z +=dz){
		norm +=pow(zdistr(z),2.0)*dz;
		array[2] = z;
		a = 1./(1+z);
		if (zphot<0.5){
			fK     = f_K(chi(a));
			k      = array[0]/fK;
			D_a=growfac(a,1.,1.)/growfac(1.,1.,1.);
			hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
			
			res= a*a*hoverh0/fK/fK;
			
			res = res*P_delta(a, k)*pow(zdistr(z),2.0);
		}
		else{res = int_GSL_integrate_qag2(int_for_C_gg_photoz,(void*)array,0.35,0.999,NULL,1000)*pow(zdistr(z),2.0);}
		RES += res*dz;
		
	}
	return cosmology.bias*cosmology.bias*RES/(norm*norm*10.0);
}

double int_for_C_gi_photoz (double a, void *params){
	double res,k,fK,D_a,GIz,hoverh0;
	double *array = (double*)params;
	
	fK     = f_K(chi(a));
	k      = array[0]/fK;
	D_a=growfac(a,1.,1.)/growfac(1.,1.,1.);
	hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
	
	res= a*a*hoverh0/fK/fK;
	
	res = res*0.0134*cosmology.Omega_m/D_a*P_delta(a, k);
	if (array[1] > 0.5){res = res*pow(p_zspec(1./a-1.,array[3])/(a*a),2.0);}
	return res;
}

double C_gi_photoz(double l, double zm, double DZ, double zphot){
	double RES,res,k,fK,D_a,GIz,hoverh0,z,a,norm,dz;
	double array[3] = {l,zphot,zm};
	dz = 0.02;
	RES = 0;
	norm = 0.;
	for (z = zm-DZ/2.0+dz/2.0; z< zm+DZ/2.0; z +=dz){
		norm +=pow(zdistr(z),2.0)*dz;
		array[2] = z;
		a = 1./(1+z);
		if (zphot<0.5){
			fK     = f_K(chi(a));
			k      = array[0]/fK;
			D_a=growfac(a,1.,1.)/growfac(1.,1.,1.);
			hoverh0 = sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
			
			res= a*a*hoverh0/fK/fK;
			
			res = res*0.0134*cosmology.Omega_m/D_a*P_delta(a, k)*pow(zdistr(z),2.0);
		}
		else{res = int_GSL_integrate_qag2(int_for_C_gi_photoz,(void*)array,0.35,0.999,NULL,1000)*pow(zdistr(z),2.0);}
		RES += res*dz;
		
	}
	return cosmology.bias*RES/(norm*norm*10.0);
}
double get_ngal (double zm, double DZ){
	double norm, z,dz;
	dz = 0.02;
	norm =0.;
	for (z = zm+DZ/2.0+dz/2.0;z<1.6; z +=dz){
		norm +=pow(zdistr(z),2.0)*dz;
	}
	return survey.n_gal*norm;
}


double get_nlens (double zm, double DZ){
	double norm,z, dz;
	dz = 0.02;
	norm =0.;
	for (z = zm-DZ/2.0+dz/2.0; z< zm+DZ/2.0; z +=dz){
		norm +=pow(zdistr(z),2.0)*dz;
	}
	return survey.n_gal*norm;
}
int main(void){
	gsl_set_error_handler_off ();
	set_survey_parameters_to_DES();
	set_redshift_BCC();
	redshift.shear_histogram_zbins=0;
	redshift.clustering_histogram_zbins=0;
	redshift.clustering_zdistrpar_zmin = 0.1;
	set_cosmological_parameters_to_WMAP_7years_BAO_SN();
	cosmology.bias=1.9;
	double k,z,a,b;
	for (z = 0.2; z< 1.2; z +=0.01){
		a = p_zphot(z,0.5);
		b = p_zspec(z,0.5);
		printf ("%e %e %e    %e %e\n",z, pf_photoz(z,-1), zdistr_photoz(z,-1), p_zspec(z,0.5), p_zphot(z,0.5));
//		printf("%e %e %e\n",k, C_gi_photoz(k,0.5,0.2,0.0),C_gi_photoz(k,0.5,0.2,1.0) );
	}
	return 0;
}

