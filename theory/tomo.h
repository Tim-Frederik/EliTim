#ifndef __TOMO_H
#define __TOMO_H

double indexcalc(int a,int b);
double zdistr_tomo(double z,int i);
double int_for_g_tomo(double aprime, void *params);
double g_source_tomo(double a, int zbin);

double int_for_p_shear_shear_tomo(double a, void *params);
double P_shear_shear_tomo(double s, int Npowerspec);
double xi_shear_shear_tomo(int pm, double theta,int Npowerspec);
void xi_via_hankel_shear_shear_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec);

double xi_shear_magnification_tomo(int pm, double theta,int Npowerspec);
void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec);

double int_for_p_shear_position_tomo(double a, void *params);
double P_shear_position_tomo(double s, int Npowerspec);
double xi_shear_position_tomo(int pm, double theta,int Npowerspec);
void xi_via_hankel_shear_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec);

double xi_magnification_position_tomo(int pm, double theta,int Npowerspec);
void xi_via_hankel_magnification_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec);

double int_for_p_position_position_tomo(double a, void *params);
double P_position_position_tomo(double s, int Npowerspec);
double xi_position_position_tomo(int pm, double theta,int Npowerspec);
void xi_via_hankel_position_position_tomo(double **xi, double *logthetamin, double *logthetamax,int Npowerspec);

#endif





























