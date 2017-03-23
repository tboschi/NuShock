#ifndef MCTOOLS_H_
#define MCTOOLS_H_

double find_emax(bool, double,double,double,double,int); 
double which_flux(bool,double, double, double, double, double, int);
double get_event(bool probornot, double mS, double Ue, double Um, double Ut, int detector_flag, gsl_rng *r, double emax);



#endif

