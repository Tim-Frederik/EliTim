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

#include "power_covariances.h"
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
extern tomopara tomo;
extern configpara configuration;
extern redshiftpara redshift;
extern globalpara global;
extern nuisancepara nuisance;
extern sur survey;




/*=========================================================*/
/*================= Matrix/Inversion routines =============*/
/*=========================================================*/


void matrix_multi(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c,int N_intern)
{
  int i,j,k;
  double d;
  for (i=0;i<N_intern;i++){
    for (j=0;j<N_intern;j++){
      d=0;
      for (k=0;k<N_intern;k++){
	d=d+gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,d);
    }
  }
}



void SVD_inversion(gsl_matrix *cov, gsl_matrix *inverseSVD,int Nmatrix)
{
  int i,j;
  gsl_matrix *V=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_matrix *U=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_vector *S=gsl_vector_calloc(Nmatrix);
  gsl_vector *work=gsl_vector_calloc(Nmatrix);
  
  gsl_matrix_memcpy(U,cov);
  
  gsl_linalg_SV_decomp(U,V,S,work);
  
  for (i=0;i<Nmatrix;i++){
    gsl_vector *b=gsl_vector_calloc(Nmatrix);
    gsl_vector *x=gsl_vector_calloc(Nmatrix);
    gsl_vector_set(b,i,1.0);
    gsl_linalg_SV_solve(U,V,S,b,x);
    for (j=0;j<Nmatrix;j++){
      gsl_matrix_set(inverseSVD,i,j,gsl_vector_get(x,j));
    }
  }
}


void eigen_sing(gsl_matrix *cov,gsl_vector *eval,gsl_vector *sing, int N_intern)
{
  gsl_matrix *mat=gsl_matrix_calloc(N_intern,N_intern);
  gsl_vector *work=gsl_vector_calloc(N_intern);
  gsl_matrix_memcpy(mat,cov);
  
  gsl_matrix *evec=gsl_matrix_calloc(N_intern,N_intern);
  gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc(2*N_intern);
  gsl_eigen_symmv(mat,eval,evec,w);
  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
  gsl_eigen_symmv_free(w);
  gsl_matrix_free(evec);
  
  gsl_matrix *mat2=gsl_matrix_calloc(N_intern,N_intern);
  gsl_matrix_memcpy(mat2,cov);
  gsl_matrix *V=gsl_matrix_calloc(N_intern,N_intern);
  gsl_linalg_SV_decomp(mat2,V,sing,work);
}


void pseudo_inverse(gsl_matrix *cov, gsl_matrix *invcov,int zero,int Nmatrix)
{
  int i,j;
  gsl_matrix *V=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_matrix *m2=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_vector *S=gsl_vector_calloc(Nmatrix);
  gsl_vector *work=gsl_vector_calloc(Nmatrix);
  gsl_matrix *sing=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_matrix *mat1=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_matrix *unity=gsl_matrix_calloc(Nmatrix,Nmatrix);
  
  gsl_matrix_memcpy(m2,cov);
  gsl_linalg_SV_decomp(m2,V,S,work);
  for (i=0;i<Nmatrix;i++){
    gsl_matrix_set(sing,i,i,1./gsl_vector_get(S,i));
    //if (i>(Nmatrix-zero)){
      //	gsl_matrix_set(sing,i,i,0.0);
    //}
    printf("pseudoinv%le\n",gsl_matrix_get(sing,i,i));
  }
  gsl_matrix_transpose(m2);
  matrix_multi(V,sing,mat1,Nmatrix);
  matrix_multi(mat1,m2,invcov,Nmatrix);
  matrix_multi(invcov,cov,unity,Nmatrix);
  for (i=0;i<Nmatrix;i++){
    for (j=0;j<Nmatrix;j++){
      if((i!=j) && (gsl_matrix_get(unity,i,j)>1.e-9)){
	printf("error in unity %d %d %le\n",i,j,gsl_matrix_get(unity,i,j));
      }
      //printf("%d %d %le\n",i,j,gsl_matrix_get(unity,i,j));
    }
  }
}


void read_write_cov_NonGauss(char *INFILE,int Ncl)
{
  double space,space1,ell1,ell2;
  int Nmatrix,k,i,j,l,m,p,q,intspace;
  FILE *F,*F2;
  char filename[600];
  
  sprintf(filename,"%s",INFILE);
  printf("%s\n",filename);
  F=fopen(filename,"r");
  sprintf(filename,"%s_easy_read",INFILE);
  printf("%s\n",filename);
  F2=fopen(filename,"w");
  Nmatrix=Ncl*tomo.shear_Npowerspectra;
  gsl_matrix *cov=gsl_matrix_calloc(Nmatrix,Nmatrix);
  k=0;
  for (i=0;i<tomo.shear_Npowerspectra; i++){
    for (j=i;j<tomo.shear_Npowerspectra; j++){
      for (p=0;p<Ncl; p++){
	for (q=p;q<Ncl; q++,k++){
	  fscanf(F,"%d %d %d %d %le %le %le %le %le %le %le %le\n",&intspace,&intspace,&intspace,&intspace,&ell1,&ell2,&space,&space,&space,&space,&space,&space1);
	  gsl_matrix_set(cov,i*Ncl+p,j*Ncl+q,space1);
	  gsl_matrix_set(cov,i*Ncl+q,j*Ncl+p,space1);
	  
	  gsl_matrix_set(cov,j*Ncl+p,i*Ncl+q,space1);
	  gsl_matrix_set(cov,j*Ncl+q,i*Ncl+p,space1);
	}
      }
    }
  }
  printf("%d\n",k);
  fclose(F);
  for (l=0;l<tomo.shear_Npowerspectra;l++){
    for (m=0;m<tomo.shear_Npowerspectra;m++){
      for (p=0;p<Ncl; p++){
	for (q=0;q<Ncl; q++){
	  fprintf(F2,"%d %d %d %d %le\n",l,m,p,q,gsl_matrix_get(cov,l*Ncl+p,m*Ncl+q));
	}
      }
    }
  }
  fclose(F2);
}

void read_write_cov_Gauss(char *INFILE, int Ncl)
{
  double space1,ell1,ell2;
  int Nmatrix,i,j,l,m,p,q,intspace;
  FILE *F,*F2;
  char filename[600];
  
  sprintf(filename,"%s",INFILE);
  printf("%s\n",filename);
  F=fopen(filename,"r");
  sprintf(filename,"%s_easy_read",INFILE);
  printf("%s\n",filename);
  F2=fopen(filename,"w");
  Nmatrix=Ncl*tomo.shear_Npowerspectra;
  gsl_matrix *cov=gsl_matrix_calloc(Nmatrix,Nmatrix);
  
  for (i=0;i<tomo.shear_Npowerspectra; i++){
    for (j=i;j<tomo.shear_Npowerspectra; j++){
      for (p=0;p<Ncl; p++){
	for (q=p;q<Ncl; q++){
	  space1=0.0;
	  if (p==q){
	    fscanf(F,"%d %d %d %d %le %le %le \n",&intspace,&intspace,&intspace,&intspace,&ell1,&ell2,&space1);
	    // printf("%le\n",space1);
	  }
	  gsl_matrix_set(cov,i*Ncl+p,j*Ncl+q,space1);
	  gsl_matrix_set(cov,i*Ncl+q,j*Ncl+p,space1);
	  
	  gsl_matrix_set(cov,j*Ncl+p,i*Ncl+q,space1);
	  gsl_matrix_set(cov,j*Ncl+q,i*Ncl+p,space1);	  
	}
      }
    }
  }
  fclose(F);
  
  for (l=0;l<tomo.shear_Npowerspectra;l++){
    for (m=0;m<tomo.shear_Npowerspectra;m++){
      for (p=0;p<Ncl; p++){
	for (q=0;q<Ncl; q++){
	  fprintf(F2,"%d %d %d %d %le\n",l,m,p,q,gsl_matrix_get(cov,l*Ncl+p,m*Ncl+q));
	}
      }
    }
  }
  fclose(F2);
}


void invcov_CL_tomo(char *INFILE, int cut_min, int cut_max,int Ncl) 
{    
  double space1;
  int Nmatrix,i,j,l,m,p,q,intspace1,intspace2,intspace3,intspace4;
  FILE *F,*F2;
  char filename[600];
  int new1=0; 
  int new2=0;
  sprintf(filename,"%s_easy_read",INFILE);
  //printf("%s\n",filename);
  F=fopen(filename,"r");
  int reduce=cut_max+cut_min;
  Nmatrix=(Ncl-reduce)*tomo.shear_Npowerspectra;
  gsl_matrix *cov=gsl_matrix_calloc(Nmatrix,Nmatrix);
  
  for (i=0;i<tomo.shear_Npowerspectra; i++){
    for (p=0;p<tomo.shear_Npowerspectra; p++){
      for (j=0;j<Ncl; j++){
	for (q=0;q<Ncl; q++){
	  fscanf(F,"%d %d %d %d %le\n",&intspace1,&intspace2,&intspace3,&intspace4,&space1);
	  if((j<Ncl-cut_max) && (q<Ncl-cut_max) && (j>=cut_min) && (q>=cut_min)){ 
	    //printf("%d %d %d %d %le %d %d\n",intspace1,intspace2,intspace3,intspace4,space1,i*(Ncl-reduce)+j-cut_min,p*(Ncl-reduce)+q-cut_min);
	    gsl_matrix_set(cov,i*(Ncl-reduce)+j-cut_min,p*(Ncl-reduce)+q-cut_min,space1);
	  }
	}
      }
    }
  }
  gsl_matrix *cop=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_vector *eigen=gsl_vector_calloc(Nmatrix);
  gsl_vector *sing=gsl_vector_calloc(Nmatrix);
  gsl_matrix_memcpy(cop,cov);
  eigen_sing(cop,eigen,sing,Nmatrix);
  
  for (i=0;i<tomo.shear_Npowerspectra*(Ncl-reduce); i++){
    //printf("%d %le %le\n",i,gsl_vector_get(eigen,i),gsl_vector_get(sing,i));
    if(gsl_vector_get(eigen,i)<0.0) {
      printf("COV %dth EV=%le\n",i,gsl_vector_get(eigen,i));
      new1=1;
    }
    //if((abs(gsl_vector_get(eigen,i)/gsl_vector_get(sing,i))-1.) >.1) printf("COV Deviation EV/SV-1 = %le\n",i,(gsl_vector_get(eigen,i)/gsl_vector_get(sing,i)-1.));
  }  
  gsl_matrix *inverseSVD=gsl_matrix_calloc(Nmatrix,Nmatrix);
  gsl_vector *eigen2=gsl_vector_calloc(Nmatrix);
  gsl_vector *sing2=gsl_vector_calloc(Nmatrix);
  
  SVD_inversion(cov,inverseSVD,Nmatrix);	
  gsl_matrix_memcpy(cop,inverseSVD);
  eigen_sing(cop,eigen2,sing2,Nmatrix);
  for (i=0;i<tomo.shear_Npowerspectra*(Ncl-reduce); i++){
    //printf("%d %le %le\n",i,gsl_vector_get(eigen,i),gsl_vector_get(sing,i));
    if(gsl_vector_get(eigen2,i)<0.0){
      printf("INV %dth EV=%le\n",i,gsl_vector_get(eigen,i));
      new2=1;
    }
    //if((gsl_vector_get(eigen2,i)/gsl_vector_get(sing2,i)-1.) >.1) printf("INV Deviation EV/SV-1 =  %le\n",i,(gsl_vector_get(eigen,i)/gsl_vector_get(sing,i)-1.));
  } 
  sprintf(filename,"%s_inverse_binmin=%d_binmax=%d",INFILE,cut_min,Ncl-cut_max);
  F2=fopen(filename,"w");
  for (l=0;l<Nmatrix;l++){
    for (m=0;m<Nmatrix;m++){
      fprintf(F2,"%le\n",gsl_matrix_get(inverseSVD,m,l));
    }
  }
  fclose(F2);
  
  if(new1==0) printf("CONGRATS! Cov Pos Def\n");
  if(new2==0) printf("CONGRATS! InvCov Pos Def\n");
  printf("condition number COV: %le\n",gsl_vector_get(sing,0)/gsl_vector_get(sing,Nmatrix-1));
  printf("condition number INV: %le\n",gsl_vector_get(sing2,0)/gsl_vector_get(sing2,Nmatrix-1));
  gsl_matrix_free(cov);
  gsl_matrix_free(inverseSVD);
  gsl_matrix_free(cop);
}



/*=========================================================*/
/*================= CovCalc routines =============*/
/*=========================================================*/



double cov_quad(double pspec[3][3], int a,int b,int c, int d)
{ 
  return pspec[a][c]*pspec[b][d]+pspec[a][d]*pspec[b][c];
}

double cov_lin(double pspec[3][3],double noise[3][3], int a,int b,int c, int d)
{ 
  return pspec[a][c]*noise[b][d]+pspec[b][d]*noise[a][c]+pspec[a][d]*noise[b][c]+pspec[b][c]*noise[a][d];
}

double cov_noise(double noise[3][3], int a,int b,int c, int d)
{ 
  return noise[a][c]*noise[b][d]+noise[a][d]*noise[b][c];
}


void make_cov_power_no_tomo(char *OUT_FILE,char *PATH,char *PARA_FILE,int N_para,int START,int END,int set_cosmo, double *ell,double *dell,int Ncl, int *calc, int NonGauss)
{
  char name[21][200]={"ssss","ssmm","sspp","sssm","sssp","ssmp","mmmm","mmpp","mmsm","mmsp","mmmp","pppp","ppsm","ppsp","ppmp","smsm","smsp","smmp","spsp","spmp","mpmp"};
  char filename[21][200];
  double prefac,A,n,noise_ss,noise_mm,noise_pp,a,b,c,fsky,prefac3,prefac2,tl[5],th1[5],th2[5],thsv[5],pre[5];
  double d,e,f,g,tot;
  int i,j,k,l; 
  double *Om,*Ob,*ns,*s8,*w0,*wa;
  FILE *ein,*F[21];
  
  Om= sm2_vector(0, N_para-1);
  Ob= sm2_vector(0, N_para-1);
  ns= sm2_vector(0, N_para-1);
  s8 = sm2_vector(0, N_para-1);
  w0= sm2_vector(0, N_para-1);
  wa= sm2_vector(0, N_para-1);
  
  ein=fopen(PARA_FILE,"r");
  for (i=0;i<N_para;i++){
    fscanf(ein,"%le %le %le %le %le %le\n",&Om[i],&Ob[i],&ns[i],&s8[i],&w0[i],&wa[i]);
  }
  fclose(ein);
  
  A=survey.area*survey.area_conversion_factor;
  n=survey.n_gal*survey.n_gal_conversion_factor;
  noise_ss= survey.sigma_e*survey.sigma_e/2.0/n;
  noise_mm= survey.sigma_e*survey.sigma_e/n;
  noise_pp= 1.0/n;
  
  fsky = survey.area/41253.0;
  
  prefac=2.0*constants.pi/A;
  prefac3= 3.0/2.0*cosmology.Omega_m;
  
  printf("surveyarea=%le sigma_e=%le n_gal=%le\n",survey.area,survey.sigma_e,survey.n_gal);
  printf("A=%le sigma_e_sqr/2/n=%le\n",A,noise_ss);
  
  double indices[21][4]={
    {0,0,0,0},{1,1,1,1},{2,2,2,2},{0,1,0,1},{0,2,0,2},{1,2,1,2},
    {0,0,1,1},{0,0,2,2},{0,0,0,1},{0,0,0,2},{0,0,1,2},
    {1,1,2,2},{1,1,0,1},{1,1,0,2},{1,1,1,2},
    {2,2,0,1},{2,2,0,2},{2,2,1,2},
    {0,1,0,2},{0,1,1,2},
    {0,2,1,2}
  };
  
  double indices_NG[21][2]={
    {4,0},{4,2},{2,0},{4,1},{3,0},{3,1},
    {4,4},{2,2},{4,3},{3,2},{3,3},
    {0,0},{2,1},{1,0},{2,1},
    {4,2},{3,1},{3,2},
    {2,0},{2,1},
    {2,2},
  };
  
  //calculating pre-factors for projected trispectra
  for(i=0;i<5;i++){
    pre[i]=pow(prefac3,1.0*i)/(2.0*constants.twopi*fsky); 
  }
  
  double table_noise[3][3],table[3][3],term_noise[21],term_lin[21],term_quad[21],c0[21];
  
  table_noise[0][0]=noise_ss;
  table_noise[0][1]=0.0;
  table_noise[0][2]=0.0;
  table_noise[1][0]=0.0;
  table_noise[1][1]=noise_mm;
  table_noise[1][2]=0.0;
  table_noise[2][0]=0.0;
  table_noise[2][1]=0.0;
  table_noise[2][2]=noise_pp;
  
  for(j=0; j<21; j++){
    c0[j]=cov_noise(table_noise,indices[j][0],indices[j][1],indices[j][2],indices[j][3]);
  }
  
  for (k=START; k<END;k++){
    printf("Om=%le Ob=%le ns=%le s8=%le w0=%le wa=%le\n",Om[k],Ob[k],ns[k],s8[k],w0[k],wa[k]);
    set_cosmology(set_cosmo);
    cosmology.Omega_m=Om[k];
    cosmology.omb=Ob[k];
    cosmology.n_spec=ns[k];
    cosmology.sigma_8=s8[k];
    cosmology.w0=w0[k];
    cosmology.wa=wa[k];
    for(l=0; l<21; l++){
      if(calc[l]==1){
	if(NonGauss==0){
	  sprintf(filename[l],"%s%sc%s_power_Gauss_no_tomo_Start%d_End%d",PATH,OUT_FILE,name[l],START,END);
	  F[l]=fopen(filename[l],"w");
	} 
	if(NonGauss==1){
	  sprintf(filename[l],"%s%sc%s_power_NonGauss_no_tomo_Start%d_End%d",PATH,OUT_FILE,name[l],START,END);
	  F[l]=fopen(filename[l],"w");
	}
      }
    }
    
    for(i=0; i<Ncl; i++){
      prefac2=prefac/ell[i]/dell[i];
      
      table[0][0]=P_shear_shear(ell[i]);
      table[0][1]=2.0*P_shear_shear(ell[i]);
      table[0][2]=P_shear_position(ell[i]);
      table[1][0]=table[0][1];
      table[1][1]=4.0*P_shear_shear(ell[i]);
      table[1][2]=2.0*P_shear_position(ell[i]);
      table[2][0]=table[0][2];
      table[2][1]=table[1][2];
      table[2][2]=P_position_position(ell[i]);
      
      for(l=0; l<21; l++){
	if(calc[l]==1){
	  term_noise[l]=prefac2*c0[l];
	  term_lin[l]=prefac2*cov_lin(table,table_noise,indices[l][0],indices[l][1],indices[l][2],indices[l][3]);
	  term_quad[l]=prefac2*cov_quad(table,indices[l][0],indices[l][1],indices[l][2],indices[l][3]);
	  if(NonGauss==0){
	    fprintf(F[l],"%le %le %le %le %le %le %le\n",ell[i]-dell[i],ell[i],ell[i]+dell[i],term_noise[l],term_lin[l],term_quad[l],term_noise[l]+term_lin[l]+term_quad[l]);
	  }
	}
      }
      if(NonGauss==1){ 
	for(j=i; j<Ncl; j++){
	  a=b=c=0.0;
	  if(i==j){
	    a=term_noise[l];
	    b=term_lin[l];
	    c=term_quad[l];
	  }
	  if (calc[0]==1){
	    tl[4]=project_tri_lin_cov(ell[i],ell[j],4); //cssss
	    th1[4]=project_tri_1h_cov(ell[i],ell[j],4);
	    th2[4]=project_tri_2h_cov(ell[i],ell[j],4);
	    thsv[4]=HSV_cov(ell[i],ell[j],fsky);
	  }
	  if (calc[2]==1){
	    tl[2]=project_tri_lin_cov(ell[i],ell[j],2); //csspp
	    th1[2]=project_tri_1h_cov(ell[i],ell[j],2);
	    th2[2]=project_tri_2h_cov(ell[i],ell[j],2);
	    thsv[2]=HSV_cov(ell[i],ell[j],fsky);
	  }
	  if (calc[4]==1){
	    tl[3]=project_tri_lin_cov(ell[i],ell[j],3); //csssp
	    th1[3]=project_tri_1h_cov(ell[i],ell[j],3);
	    th2[3]=project_tri_2h_cov(ell[i],ell[j],3);
	    thsv[3]=HSV_cov(ell[i],ell[j],fsky);
	  }
	  if (calc[11]==1){
	    tl[0]=project_tri_lin_cov(ell[i],ell[j],0); //cpppp
	    th1[0]=project_tri_1h_cov(ell[i],ell[j],0);
	    th2[0]=project_tri_2h_cov(ell[i],ell[j],0);
	    thsv[0]=HSV_cov(ell[i],ell[j],fsky);
	  }
	  if (calc[13]==1){
	    tl[1]=project_tri_lin_cov(ell[i],ell[j],0); //csppp
	    th1[1]=project_tri_1h_cov(ell[i],ell[j],0);
	    th2[1]=project_tri_2h_cov(ell[i],ell[j],0);
	    thsv[1]=HSV_cov(ell[i],ell[j],fsky);
	  }
	  for(l=0; l<21; l++){
	    if (calc[l]==1){
	      int alpha=indices_NG[l][1];
	      int beta=indices_NG[l][2];
	      d=pow(2.0,beta)*pre[alpha]*tl[alpha];
	      e=pow(2.0,beta)*pre[alpha]*th1[alpha];
	      f=pow(2.0,beta)*pre[alpha]*th2[alpha];
	      g=thsv[alpha];
	      tot=a+b+c+d+e+f+g;
	      fprintf(F[l],"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",ell[i]-dell[i],ell[i],ell[i]+dell[i],ell[j]-dell[j],ell[j],ell[j]+dell[j],a,b,c,d,e,f,g,tot);
	      
	    }
	  }
	}
      }
    }
    for(l=0; l<21; l++){
      if (calc[l]==1){
	fclose(F[l]);
      }	
    }
  }
}



void make_cov_power_tomo(char *OUT_FILE,char *PATH, char *PARA_FILE,int N_para,int START,int END,int set_cosmo, double *ell,double *dell,int Ncl,int *calc, int NonGauss)
{
  char filename[400];
  double prefac,A,ngal,noise_ss,fsky,prefac3,prefac2,tl,th1,th2,thsv,a,gauss;
  double array[6];
  int i,j,k,l,m,n,o,i1,i2,i3,i4;
  FILE *F1,*ein;
  
  double Pi1,Pi2,Pi3,Pi4;
  double *Om,*Ob,*ns,*s8,*w0,*wa;
  
  Om= sm2_vector(0, N_para-1);
  Ob= sm2_vector(0, N_para-1);
  ns= sm2_vector(0, N_para-1);
  s8 = sm2_vector(0, N_para-1);
  w0= sm2_vector(0, N_para-1);
  wa= sm2_vector(0, N_para-1);
  
  ein=fopen(PARA_FILE,"r");
  for (i=0;i<N_para;i++){
    fscanf(ein,"%le %le %le %le %le %le\n",&Om[i],&Ob[i],&ns[i],&s8[i],&w0[i],&wa[i]);
  }
  fclose(ein);
  
  A=survey.area*survey.area_conversion_factor;
  ngal=survey.n_gal*survey.n_gal_conversion_factor;
  ngal=ngal/tomo.shear_Nbin; 
  noise_ss= survey.sigma_e*survey.sigma_e/2.0/ngal;
  
  fsky = survey.area/41253.0;
  
  prefac=2.0*constants.pi/A;
  prefac3= 3.0/2.0*cosmology.Omega_m;
  printf("surveyarea=%le sigma_e=%le n_gal=%le\n",survey.area,survey.sigma_e,survey.n_gal);
  printf("A=%le sigma_e_sqr/2/n=%le\n",A,noise_ss);
  
  
  for (k=START; k<END;k++){
    printf("Om=%le Ob=%le ns=%le s8=%le w0=%le wa=%le\n",Om[k],Ob[k],ns[k],s8[k],w0[k],wa[k]);
    set_cosmology(set_cosmo);
    cosmology.Omega_m=Om[k];
    cosmology.omb=Ob[k];
    cosmology.n_spec=ns[k];
    cosmology.sigma_8=s8[k];	
    cosmology.w0=w0[k];
    cosmology.wa=wa[k];
    if(NonGauss==0){
      sprintf(filename,"%s%s_power_lmin=%.3le_lmax=%.3le_Ncl=%d_tomo_Gauss",PATH,OUT_FILE,ell[0],ell[Ncl-1],Ncl);
      F1=fopen(filename,"w");
    }
    
    if(NonGauss==1){
      sprintf(filename,"%s%s_power_lmin=%.3le_lmax=%.3le_Ncl=%d_tomo_NonGauss",PATH,OUT_FILE,ell[0],ell[Ncl-1],Ncl);
      F1=fopen(filename,"w");
    }
    
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=l;m<tomo.shear_Nbin;m++){
	for (n=l;n<tomo.shear_Nbin; n++){
	  for (o=n;o<tomo.shear_Nbin; o++){
	    if((o>=m)|| (n>l) ){  //both important 1)otherwise lower triangular matrix is sampled 2) otherwise e.g. 15 23 is excluded
	     i1=indexcalc(l,n);
	     i2=indexcalc(m,o);
	     i3=indexcalc(l,o);
	     i4=indexcalc(m,n);
	     for(i=0; i<Ncl; i++){
	       Pi1=P_shear_shear_tomo(ell[i],i1);
	       Pi2=P_shear_shear_tomo(ell[i],i2); 
	       Pi3=P_shear_shear_tomo(ell[i],i3); 
	       Pi4=P_shear_shear_tomo(ell[i],i4); 
	       //printf("%le %le %le %le\n",Pi1,ell[i],dell[i],prefac);
	       if(l==n) Pi1=Pi1+noise_ss;
	       if(m==o) Pi2=Pi2+noise_ss;
	       if(l==o) Pi3=Pi3+noise_ss;
	       if(m==n) Pi4=Pi4+noise_ss;
	       
	       prefac2=prefac/ell[i]/dell[i];
	       gauss=prefac2*(Pi1*Pi2+Pi3*Pi4);
	       if(NonGauss==0){
		 fprintf(F1,"%d %d %d %d %le %le %le\n",l,m,n,o,ell[i],ell[i],gauss);
	       }
	       if(NonGauss==1){
		 for(j=i; j<Ncl; j++){
		   a=0.0;
		   if(i==j){
		     a=gauss;
		   }
		   //printf("%le %le %le\n",a,b,c);
		   array[0] = ell[i];
		   array[1] = ell[j];
		   array[2] = (double) l;
		   array[3] = (double) m;
		   array[4] = (double) n;
		   array[5] = (double) o;
		   tl = pow(prefac3,4.0)*int_GSL_integrate_crude(inner_project_tri_lin_cov_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax +1.),0.999999,NULL,1000)/(2.0*constants.twopi*fsky);
		   //printf("tl %le\n",tl);
		   th1 = pow(prefac3,4.0)*int_GSL_integrate_crude(inner_project_tri_1h_cov_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax +1.),0.999999,NULL,1000)/(2.0*constants.twopi*fsky);
		   //printf("th1 %le\n",th1);
		   th2 = pow(prefac3,4.0)*int_GSL_integrate_crude(inner_project_tri_2h_cov_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax +1.),0.999999,NULL,1000)/(2.0*constants.twopi*fsky);
		   //printf("th2 %le\n",th2);
		   thsv=HSV_cov_tomo(ell[i],ell[j],fsky,l,m,n,o);
		   fprintf(F1,"%d %d %d %d %le %le %le %le %le %le %le %le\n",l+1,m+1,n+1,o+1,ell[i],ell[j],a,tl,th1,th2,thsv,a+tl+th1+th2+thsv);
		   printf("%d %d %d %d %le %le %le %le %le %le %le %le\n",l+1,m+1,n+1,o+1,ell[i],ell[j],a,tl,th1,th2,thsv,a+tl+th1+th2+thsv);
		 }
	       }
	     }
	    }
	  }
	}
      }
    }
    fclose(F1);
    if (NonGauss==1)read_write_cov_NonGauss(filename,Ncl);
    if(NonGauss==0) read_write_cov_Gauss(filename,Ncl);
    invcov_CL_tomo(filename,0,0,Ncl); 
  }  
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


void write_parameters(char *PARA_FILE, char*PATH, char *OUT_FILE,int N_para,int START,int END,int shear_shear,int mag_mag,int pos_pos,int shear_mag,int shear_pos,int mag_pos,int Ncl, double lmax, double lmin,int NonGauss)
{
  printf("---------------------------------------\n");
  printf("---------RUN MODE PARAMETERS-----------\n");
  printf("\n");
  printf("PARA_FILE: %s\n",PARA_FILE);
  printf("PATH: %s\n",PATH);
  printf("OUT_FILE: %s\n",OUT_FILE);
  printf("TOTAL # cosmologies: %d\n",N_para);
  printf("START: %d\n",START);
  printf("END: %d\n",END);
  printf("\n");
  
  if (shear_shear==1) printf("# Computing shear shear power spectrum covariance\n");
  if (mag_mag==1) printf("# Computing mag mag power spectrum covariance\n");
  if (pos_pos==1) printf("# Computing pos pos power spectrum covariance\n");
  if (shear_mag==1) printf("# Computing shear mag power spectrum covariance\n");
  if (shear_pos==1) printf("# Computing shear pos power spectrum covariance\n");
  if (mag_pos==1) printf("# Computing mag pos power spectrum covariance\n");
  
  if(NonGauss==1)     printf("NonGauss: YES\n");
  if(NonGauss==0)     printf("NonGauss: NO\n");
  
  printf("# bins CL: %d\n",Ncl);
  printf("lmax: %le\n",lmax);
  printf("lmin: %le\n",lmin);
  
  printf("---------------------------------------\n");
  printf("---------REDSHIFT PARAMETERS-----------\n");
  printf("\n");
  
  if(configuration.TOMOGRAPHY==0) {
    printf("Tomography: NO\n");
    printf("z_min: %le\n",redshift.shear_zdistrpar_zmin);
    printf("z_max: %le\n",redshift.shear_zdistrpar_zmax);
  }
  if(configuration.TOMOGRAPHY==1) {
    printf("Tomography: YES %d z-bins %d power spectra\n",tomo.shear_Nbin,tomo.shear_Npowerspectra);
    printf("z1: %le - %le\n",tomo.shear_zmin[0],tomo.shear_zmax[0]);
    printf("z2: %le - %le\n",tomo.shear_zmin[1],tomo.shear_zmax[1]);
    printf("z3: %le - %le\n",tomo.shear_zmin[2],tomo.shear_zmax[2]);
    printf("z4: %le - %le\n",tomo.shear_zmin[3],tomo.shear_zmax[3]);
    printf("z5: %le - %le\n",tomo.shear_zmin[4],tomo.shear_zmax[4]);
  }
  
  if(redshift.shear_histogram_zbins==0){ 
    printf("Using Redshift Parametrization\n");
    printf("shear z0=%le\n",redshift.shear_z0);
    printf("shear beta: %le\n",redshift.shear_beta_p);
    printf("shear alpha: %le\n", redshift.shear_alpha);
  }
  if(redshift.shear_histogram_zbins!=0){ 
    printf("USING SHEAR REDSHIFT_FILE: %s\n",redshift.shear_REDSHIFT_FILE);
    printf("N_bins in dn/dz: %d\n",redshift.shear_histogram_zbins);
  }
}    


int main(int argc, char** argv)
{
  gsl_set_error_handler_off ();
  int i,j,k;
  int probe[6],calc[21];
  char PARAFILE[200],OUTFILE[200],PATH[200];
  // RUN MODE setup
  sprintf(PARAFILE,"%s",argv[1]);
  sprintf(PATH,"%s",argv[2]);
  sprintf(OUTFILE,"%s",argv[3]);
  int N_para=atoi(argv[4]);
  int START=atoi(argv[5]);
  int END=atoi(argv[6]);
  
  int shear_shear=probe[0]=atoi(argv[7]); 
  int mag_mag=probe[1]=atoi(argv[8]); 
  int pos_pos=probe[2]=atoi(argv[9]); 
  int shear_mag=probe[3]=atoi(argv[10]); 
  int shear_pos=probe[4]=atoi(argv[11]); 
  int mag_pos=probe[5]=atoi(argv[12]); 
  
  int Ncl=atoi(argv[13]);
  double lmax= atof(argv[14]);
  double lmin= atof(argv[15]);
  
  // COSMOLOGY CODE setup
  configuration.transferfunction_EH99=atoi(argv[16]);
  configuration.COYOTE_UNIVERSE_CALIBRATION=atoi(argv[17]);
  int set_cosmo=atoi(argv[18]);
  set_cosmology(set_cosmo);
  // Survey  setup
  set_survey_parameters_to_DES();
  survey.area=atof(argv[19]);
  survey.n_gal=atof(argv[20]);
  survey.sigma_e=atof(argv[21]);
  
  
  // REDSHIFT setup
  redshift.shear_photoz=atoi(argv[22]);
  redshift.shear_zdistrpar_zmin=atof(argv[23]);
  redshift.shear_zdistrpar_zmax=atof(argv[24]);
  redshift.shear_histogram_zbins=atoi(argv[25]);
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",argv[26]);
  
  redshift.shear_z0=atof(argv[27]); 
  redshift.shear_alpha=atof(argv[28]); 
  redshift.shear_beta_p=atof(argv[29]); 
  
  redshift.clustering_photoz=atoi(argv[30]);
  redshift.clustering_zdistrpar_zmin=atof(argv[31]);
  redshift.clustering_zdistrpar_zmax=atof(argv[32]);
  redshift.clustering_histogram_zbins=atoi(argv[33]);
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",argv[34]);
  
  redshift.clustering_z0=atof(argv[35]); 
  redshift.clustering_alpha=atof(argv[36]); 
  redshift.clustering_beta_p=atof(argv[37]); 
 
  
  configuration.TOMOGRAPHY=atoi(argv[38]); 
  int NonGauss=atoi(argv[39]);
  if(configuration.TOMOGRAPHY==1){
    tomo.shear_Nbin=atoi(argv[40]);
    if (tomo.shear_Nbin>10){ 
      printf("ERROR: 10 shear tomography bins MAX!");
    }
    for (i=0;i<tomo.shear_Nbin; i++){
      tomo.shear_Npowerspectra=tomo.shear_Npowerspectra+i+1;
    }
    for (i=0;i<tomo.shear_Nbin; i++){
      tomo.shear_zmin[i]=atof(argv[41+i]);
      tomo.shear_zmax[i]=atof(argv[42+i]);
      //     printf("%le %le\n",tomo.shear_zmax[i],tomo.shear_zmin[i]);
    }
    tomo.clustering_Nbin=atoi(argv[42+tomo.shear_Nbin]);
    if (tomo.clustering_Nbin>10){ 
      printf("ERROR: 10 clustering tomography bins MAX!");
    }
    for (i=0;i<tomo.clustering_Nbin; i++){
      tomo.clustering_Npowerspectra=tomo.clustering_Npowerspectra+i+1;
    }
    for (i=0;i<tomo.clustering_Nbin; i++){
      tomo.clustering_zmin[i]=atof(argv[43+tomo.shear_Nbin+i]);
      tomo.clustering_zmax[i]=atof(argv[44+tomo.shear_Nbin+i]);
      //     printf("%le %le\n",tomo.clusteringzmax[i],tomo.clusteringzmin[i]);
    }
  }
   
  write_parameters(PARAFILE,PATH,OUTFILE,N_para,START,END,shear_shear,mag_mag,pos_pos,shear_mag,shear_pos,mag_pos,Ncl,lmax,lmin,NonGauss);
  
  //binning  calculation
  //  printf("%le %le\n",lmin,lmax);
  double logdl=(log(lmax)-log(lmin))/Ncl;
  double *ell, *dell;
  ell=sm2_vector(0,Ncl-1);
  dell=sm2_vector(0,Ncl-1);
  for(i=0; i<Ncl ; i++){
    ell[i]=exp(log(lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(lmin)+(i+1)*logdl) - exp(log(lmin)+(i*logdl));
    //  printf("%le %le\n",ell[i],dell[i]);
  }
  k=0;
  for(i=0; i<6; i++){
    for(j=i; j<6; j++,k++){
      calc[k]=0;
      if(probe[i]==1 && probe[j]==1 ) calc[k]=1;
      //printf("%d\n",calc[k]);
    }    
  }
  
   //Start of routines
  if((configuration.TOMOGRAPHY==0)){ make_cov_power_no_tomo(OUTFILE,PATH,PARAFILE,N_para,START,END,set_cosmo,ell,dell,Ncl,calc,NonGauss);
  //name, min, max
  
  }
  if((configuration.TOMOGRAPHY==1)){ make_cov_power_tomo(OUTFILE,PATH,PARAFILE,N_para,START,END,set_cosmo,ell,dell,Ncl,calc,NonGauss);
  }
  
  return 0;
}


