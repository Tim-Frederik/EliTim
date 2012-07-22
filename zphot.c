#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "maths.h"
#include "zphot.h"

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

int main(){
	double z;
	double a1,a2,a3;
	for (z =0.0; z<2.0;z+=.05){
		a1 += p_zspec(z,0.2)*0.05;
		a2 += p_zspec(z,0.8)*0.05;
		a3 += p_zspec(z,1.5)*0.05;
	}
	printf("%e %e %e\n",a1,a2,a3);
	a1 = 0.; a2 = 0.; a3 = 0.;
	for (z =0.0; z<2.0;z+=.05){
		a1 += p_zphot(z,0.2)*0.05;
		a2 += p_zphot(z,0.8)*0.05;
		a3 += p_zphot(z,1.5)*0.05;
	}
	printf("%e %e %e\n",a1,a2,a3);
}