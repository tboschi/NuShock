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

#define NUMEVENTS 50000



/* ########## Main function ############### */
int main(int argc, char* argv[])
{

const gsl_rng_type * T; // Standard invocation prayer to the RNG
gsl_rng * r;
gsl_rng_env_setup();

T = gsl_rng_default;
r = gsl_rng_alloc (T);



const struct option longopts[] = 
{
	{"using-muboone", 	no_argument, 		0, '1'},
	{"using-icarus", 	no_argument, 		0, '2'},
	{"using-sbnd", 		no_argument, 		0, '3'},
	{"using-ps191",		no_argument, 		0, '4'},
	{"using-all-sbn", 	no_argument, 		0, '5'},
	{"ue4", 		no_argument, 		0, 'E'},
	{"um4", 		no_argument, 		0, 'U'},
	{"ut4", 		no_argument, 		0, 'T'},
	{"prob",		no_argument,		0, 'P'},
	{0,			no_argument, 		0,  0},
};

int index; 
int iarg = 0;
opterr=1;

bool Ue_flag = false;
bool Um_flag = false;
bool Ut_flag = false;
int detector_flag = DET_SBND;

bool decay_flag = false;

double mS = 0.0010;

while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "m:12345EUTP", longopts, &index);

	switch(iarg)
	{
		case 'm':
			mS = strtof(optarg,NULL);
			break;
		case '1':
			detector_flag = DET_MUBOONE;
			break;
		case '2':
			detector_flag = DET_ICARUS;
			break;	
		case '3':
			detector_flag = DET_SBND;
			break;
		case '4':
			detector_flag = DET_PS191;
			break;
		case '5':
			detector_flag = DET_ALL_SBN;
			break;
		case 'E':
			Ue_flag=true;	
			break;
		case 'U':
			Um_flag=true;
			break;
		case 'T':
			Ut_flag=true;
			break;
		case 'P':
			decay_flag=true;
			break;
		case '?':
			std::cout<<"Abandon hope all ye who enter this value. "<<std::endl<<std::endl;
			std::cout<<"Allowed arguments:"<<std::endl;
			std::cout<<"\t-m\t\t\tsets the sterile mass. [default = 0.0010]"<<std::endl;
			std::cout<<"\t--ue4\t\t\tToggle on/off Ue4 mixing. [default off]"<<std::endl;
			std::cout<<"\t--um4\t\t\tToggle on/off Um4 mixing. [default on]"<<std::endl;
			std::cout<<"\t--ut4\t\t\tToggle on/off Ut4 mixing. [default off]"<<std::endl;
			std::cout<<"\t--using-sbnd\t\truns for SBND."<<std::endl;
			std::cout<<"\t--using-muboone\t\truns for muBooNE."<<std::endl;
			std::cout<<"\t--using-icarus\t\truns for ICARUS."<<std::endl;
			std::cout<<"\t--using-ps191\t\truns for PS191."<<std::endl;
			std::cout<<"\t--using-all-sbn\t\truns for all three detectors of SBN."<<std::endl;
			std::cout<<"\t--prob\t\t sampels from flux/E to sample decay prob"<<std::endl;
			return 0;
	}

}

double Ue = 0.0;
double Um = 0.0;
double Ut = 0.0;

if(Ue_flag == true)
{
	Ue = 1.0;
}
else if(Um_flag == true)
{
	Um = 1.0;
}
else if(Ut_flag == true)
{
	Ut = 1.0;
}

double E=0.001; //GeV energy of particle.

double emax = find_emax(decay_flag,mS,Ue,Um,Ut,detector_flag);
//std::cout<<"emax is "<<emax<<std::endl;

static double flux_array[NUMEVENTS];
double height_draw = 0.0;
double norm_flux = 0.0;

double current = 0;




int total = 0;
while(total < NUMEVENTS)
{
	E = 5.0*gsl_rng_uniform(r);
	height_draw = gsl_rng_uniform(r);
	norm_flux = which_flux(decay_flag,E,mS,Ue,Um,Ut,detector_flag)/emax;
	if(height_draw <= norm_flux)
	{
	//	std::cout<<total<<" "<<E<<" "<<height_draw<<" "<<norm_flux<<" "<<emax<<std::endl;
		std::cout<<E<<"\t0.999999"<<std::endl;
		flux_array[total] = E;
		total++;
	}
}

//int n=0;
//#for(n=0;n<NUMEVENTS;n++)
//{
//	std::cout<<flux_array[n]<<"\t0.999999"<<std::endl;
//}

return 0;
}// end main


