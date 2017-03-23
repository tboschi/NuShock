#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define no_argument 0
#define required_argument 1
#define optional_argument 2

#include "flux.h"	  // This includes the detector-specific cut functions.
#include "mctools.h"

double get_event(bool probornot, double mS, double Ue, double Um, double Ut, int detector_flag, gsl_rng *r, double emax)
{
	double height_draw = 0.0;
	double norm_flux = 0.0;
	double E = 0.001;
	double current = 0;

	E = 5.0*gsl_rng_uniform(r);
	height_draw = gsl_rng_uniform(r);
	norm_flux = which_flux(probornot,E,mS,Ue,Um,Ut,detector_flag)/emax;
	while(height_draw > norm_flux)
	{
		E = 5.0*gsl_rng_uniform(r);
		height_draw = gsl_rng_uniform(r);
		norm_flux = which_flux(probornot, E,mS,Ue,Um,Ut,detector_flag)/emax;
	
	}

return E;
}






double find_emax(bool decay_flag, double mS, double Ue, double Um, double Ut, int detector_flag)
{

double emax = 0.0;
double current = 0.0;

for(double iE = 1e-4; iE<5+1e-5 ; iE = iE+1e-3)
{
	current = which_flux(decay_flag,iE,mS,Ue,Um,Ut,detector_flag);
//	std::cout<<iE<<" "<<current<<std::endl;
	if(current > emax)
	{
		emax = current;
//		std::cout<<iE<<" loop "<<emax<<std::endl;
	}

}


//for(E=1e-4;E<10.0+1e-5;E+=1e-3)
//{
//	current = which_flux(probornot,E,mS,Ue,Um,Ut,detector_flag);
//	std::cout<<E<<" "<<current/max<<std::endl;
//}
//	std::cout<<"internal emax is: "<<emax<<std::endl;
return emax;
}




double which_flux(bool probornot, double E, double ms, double Ue, double Uu,double Ut, int det)
{

	double ans =0.0;
       
	if(E>ms)
	{
		ans =	flux(E,ms,Ue,Uu,Ut,det);

		if(probornot)
		{
			ans= flux(E,ms,Ue,Uu,Ut,det)/sqrt(E*E-ms*ms);
		}
	}

return ans;	

}
