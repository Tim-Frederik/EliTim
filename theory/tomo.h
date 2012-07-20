#ifndef __TOMO_H
#define __TOMO_H

double zdistr_tomo1(double z);
double zdistr_tomo2(double z);
double int_for_g_tomo1(double aprime);
double int_for_g_tomo2(double aprime);
double g_source_tomo(double a, int pm);
double int_for_p_shear_shear_tomo(double a);

double P_shear_shear_tomo(double s);
double xi_shear_shear_tomo(int pm, double theta);
void xi_via_hankel_shear_shear_tomo(double **xi, double *logthetamin, double *logthetamax);
double xi_shear_magnification_tomo(int pm, double theta);
void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax);
double xi_shear_magnification_tomo(int pm, double theta);
void xi_via_hankel_shear_magnification_tomo(double **xi, double *logthetamin, double *logthetamax);

double pf_tomo1(double z);
double pf_tomo2(double z);

double int_for_p_shear_position_tomo(double a);
double P_shear_position_tomo(double s);
double xi_shear_position_tomo(int pm, double theta);
void xi_via_hankel_shear_position_tomo(double **xi, double *logthetamin, double *logthetamax);

double xi_magnification_position_tomo(int pm, double theta);
void xi_via_hankel_magnification_position_tomo(double **xi, double *logthetamin, double *logthetamax);

double int_for_p_position_position_tomo(double a);
double P_position_position_tomo(double s);
double xi_position_position_tomo(int pm, double theta);
void xi_via_hankel_position_position_tomo(double **xi, double *logthetamin, double *logthetamax);


#endif





























