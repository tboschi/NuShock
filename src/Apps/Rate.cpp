#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <getopt.h>

#include "EventGenerator.h"
#include "FluxDriver.h"
#include "DecayRates.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TGraph.h"

void Usage(char* argv0);

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"xsection", 	required_argument, 	0, 'x'},
		{"flux", 	required_argument,	0, 'f'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string DetConfig;
	TFile *XsecFile, *FluxFile;
	std::ofstream OutFile;
	
	while((iarg = getopt_long(argc,argv, "d:x:f:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'x':
				XsecFile = new TFile(optarg, "OPEN");
				break;
			case 'f':
				FluxFile = new TFile(optarg, "OPEN");
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
				return 1;
			default:
				break;
		}
	}

	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	Detector * TheBox = new Detector(DetConfig, NULL);

	TH1D* hFlux = dynamic_cast<TH1D*> (FluxFile->Get("ptotal"));
	hFlux->Scale(1.0/pow(TheBox->GetElement("Baseline"), 2));
	//hFlux->Scale(TheBox->GetElement("POT/s"));
	hFlux->Scale(TheBox->GetElement("Height")*TheBox->GetElement("Width")*1.0e4);

	TIter Next(XsecFile->GetListOfKeys());
	TKey *Kkk = dynamic_cast<TKey*> (Next());
	TDirectory *Dir = dynamic_cast<TDirectory*> (Kkk->ReadObj());
	TGraph * gXsecCC = dynamic_cast<TGraph*> (Dir->Get("tot_cc"));
	TGraph * gXsecNC = dynamic_cast<TGraph*> (Dir->Get("tot_nc"));

	TH1D * hXsec = new TH1D("hxsec", "hxsec", hFlux->GetNbinsX(), 0, 20);

	double x = 0, y = 0;
	unsigned int j = 0;
	for (unsigned int i = 1; i < hFlux->GetNbinsX()+1; ++i)
	{
		double sum = 0;
		unsigned int N = 0;
		while (x < hFlux->GetBinCenter(i) + 0.5*20.0/hFlux->GetNbinsX())
		{
			gXsecCC->GetPoint(j, x, y);
			sum += y;
			++N;
			gXsecNC->GetPoint(j, x, y);
			sum += y;
			++N;
			++j;
		}
		hXsec->SetBinContent(i, 1.0e-38*sum/N);
	}

	double MolarMass = 1;
	if (TheBox->GetElement("InTarget") == 1 || TheBox->GetElement("InTarget") == 2)	//Argon
		MolarMass = 40.0;	//grams
	else if (TheBox->GetElement("InTarget") == 3)	//Iron
		MolarMass = 56.0;	//grams

	double Ntarget = Const::fNa * TheBox->GetElement("Weight") * 1e6 / MolarMass;	//number of targets

	double Rate = 0;
	double Npot = 0;
	for (unsigned int i = 0; i < hFlux->GetNbinsX(); ++i)
	{
		Npot += hFlux->GetBinContent(i) * 0.2;
		Rate += hFlux->GetBinContent(i) * hXsec->GetBinContent(i) * 0.2;
	}
	Rate *= Ntarget;
	Rate *= TheBox->GetElement("Fiducial");		//fiducial volume

	FluxFile->Close();
	XsecFile->Close();

	double Ysec = 1.0e7 * TheBox->GetElement("Years");
	Out << "Total number of events is " << std::setprecision(10) << Rate*Ysec << " in " << TheBox->GetElement("Years") << " years" << std::endl;
	Out << "Frequency " << Rate << " Hz on average" << std::endl;
	Out << "Scale factor " << Rate*Ysec * 1e-6 << std::endl;
	Out << "tot numb of neutrinos per POT " << Npot << std::endl;

	return 0;
}
	

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -x,  --xsection" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --flux" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file, to save message" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
