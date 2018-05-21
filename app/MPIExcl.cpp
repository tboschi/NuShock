#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <getopt.h>
#include <set>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"
#include "Detector.h"

#include "TH2D.h"
#include "TFile.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument, 	0, 's'},
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"threshold", 	required_argument,	0, 't'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string SMConfig, DetConfig;
	std::string FluxConfig, EffFile;
	std::ofstream OutFile;
	//TFile *OutFile;
	std::string Channel = "ALL";
	bool UeFlag = false;
	bool UmFlag = false;
	bool UtFlag = false;
	bool Efficiency = false;
	
	while((iarg = getopt_long(argc,argv, "s:d:f:c:t:o:WEMTh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 's':
				SMConfig.assign(optarg);
				break;
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'f':
				FluxConfig.assign(optarg);
				break;
			case 'c':
				Channel.assign(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				//OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'W':
				Efficiency = true;
				break;
			case 'E':
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//MPI STUFF
	int ierr, id, np;
       	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank( MPI_COMM_WORLD, &id);
	ierr = MPI_Comm_rank( MPI_COMM_WORLD, &np);

	std::stringstream sProcId;
	sProcId << id;
	OutFile += sProcId.str();	//append proc ID to each outfile

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	EventGenerator * EvGen = new EventGenerator(SMConfig, DetConfig, FluxConfig);
	
	if (UeFlag)
		EvGen->SetChannel(Channel, Efficiency, 'E');
	if (UmFlag)
		EvGen->SetChannel(Channel, Efficiency, 'M');

	EvGen->SetMass(0);
	EvGen->SyncUu();
	
	unsigned int Grid = 500;

	double logMassStart = -2.0;
	double logMassEnd   = -0.3;
	double lmA = ( logMassStart * (np - id) + logMassEnd * id ) / np;
	double lmB = ( logMassStart * (np - id - 1) + logMassEnd * (id + 1) ) / np;
	double lmStep = (lmB - lmA) / (Grid / p);

	double logUuStart = -20.0;
	double logUuEnd   =   0.0;
	double logUuStep  = (logUuEnd - logUuStart) / Grid;

	double Mass, UuA, UuB;
	double contMass, contUu, contN;
	std::vector<double> vSignal;	//summing over energy, array of Uus
	std::vector<bool> vAbove;	//summing over energy, array of Uus
	unsigned int iMin;
	
	for (double logMass = -2.0; logMass < -0.3; logMass += 1.7/Grid)	//increase mass log
	{
		Mass = pow(10.0, logMass);
		std::cout << "Mass " << Mass << std::endl;
		EvGen->SetMass(Mass);
		EvGen->MakeFlux(1);


		for (double logUuA = -20.0; logUuA+1e-6 < 0.0; logUuA += 20.0/Grid)	//increase Uu logarithmically
		{
			UuA = pow(10.0, 0.5*logUuA);
			if (UeFlag)
				EvGen->SetUe(UuA);
			else if (UmFlag)
				EvGen->SetUm(UuA);
			//if (UtFlag)
			//	EvGen->SetUt(UuA);

			vSignal.clear();
			vSignal.resize(Grid);	//number of Uus probing
			vAbove.clear();
			vAbove.resize(Grid);	//number of Uus probing

			//iMin = Grid+1;
			iMin = 0;
			double Start, End;
			double EnStep = EvGen->GetRange(Start, End)/EvGen->GetBinNumber();
			for (double Energy = Start; Energy < End; Energy += EnStep)
			{
				unsigned int i = 0;
				for (double logUuB = -20.0; logUuB+1e-6 < 0.0; logUuB += 20.0/Grid, ++i)	//increase Uu logarithmically
				{
					UuB = pow(10.0, 0.5*logUuB);
					//if (UeFlag)
					//	EvGen->SetUe(UuB);
					if (UmFlag)
						EvGen->SetUm(UuB);
					else if (UtFlag)
						EvGen->SetUt(UuB);
					//EvGen->SetUe(Uu, 1);	//production
					//EvGen->SetUm(Uu, 0);	//decay
	
					vSignal.at(i) += EnStep * EvGen->DecayNumber(Energy, Efficiency);
				}
			}

			unsigned int j = 0;
			for (double logUuB = -20.0; logUuB+1e-6 < 0.0; logUuB += 20.0/Grid, ++j)	//increase Uu logarithmically
			{
				UuB = pow(10.0, 0.5*logUuB);
				//Out << Mass << "\t" << Uu*Uu << "\t" << vSignal.at(j) << std::endl;
				Out << Mass << "\t" << UuA*UuA << "\t" << UuB*UuB << "\t" << vSignal.at(j) << std::endl;
			}
		}
	}

	//OutFile->cd();
	//logCont->Write();
	//Contour->Write();
	//OutFile->Close();

	return 0;
}
	
void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -s,  --smconfig" << std::endl;
	std::cout << "\t\tStandard Model configuration file" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -c,  --channel" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
