#ifndef __SURVEY_H
#define __SURVEY_H

typedef struct {
     double area;		/* survey_area in deg^2. */
     double n_gal;		/* galaxy density per arcmin^2 */
     double sigma_e;		/* rms inrinsic ellipticity noise*/
     double area_conversion_factor; /*factor from deg^2 to radian^2: 60*60*constants.arcmin*constants.arcmin */
     double n_gal_conversion_factor; /*factor from n_gal/arcmin^2 to n_gal/radian^2: 1.0/constants.arcmin/constants.arcmin */
}sur;

void set_survey_parameters_to_CFHTLS();
void set_survey_parameters_to_DES();
void set_survey_parameters_to_LSST();
void set_survey_parameters_to_EUCLID();

#endif
