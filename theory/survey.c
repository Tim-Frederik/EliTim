#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "survey.h"
#include "maths.h"
extern con constants;

sur survey = {
0.0,
0.0,
0.0,
0.0,
0.0,
};

void set_survey_parameters_to_CFHTLS()
{
	survey.area   = 170.0;
	survey.n_gal   = 13.3;
	survey.sigma_e   = 0.42;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_DES()
{
	survey.area   = 5000.0;
	survey.n_gal   = 12.0;
	survey.sigma_e   = 0.3;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}


void set_survey_parameters_to_LSST()
{
	survey.area   = 15000.0;
	survey.n_gal   = 30.0;
	survey.sigma_e   = 0.22;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_EUCLID()
{
	survey.area   = 20000.0;
	survey.n_gal   = 40.0;
	survey.sigma_e   = 0.22;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_DES_DC5()
{
	survey.area   = 200.0;
	survey.n_gal   = 12.0;
	survey.sigma_e   = 0.35;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

void set_survey_parameters_to_SATO()
{
	survey.area   = 25.0;
	survey.n_gal   = 30.0;
	survey.sigma_e   = 0.22;
	survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
	survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}

