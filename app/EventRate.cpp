#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <getopt.h>

#include "flux.h"
#include "detector.h"

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
	std::string DetConfig, module;
	TFile *XsecFile, *FluxFile;
	std::ofstream OutFile;
	std::string type;
	
	while((iarg = getopt_long(argc,argv, "d:l:x:f:t:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				DetConfig.assign(optarg);
				break;
			case 'l':
				module.assign(optarg);
				break;
			case 'x':
				XsecFile = new TFile(optarg, "OPEN");
				break;
			case 'f':
				FluxFile = new TFile(optarg, "OPEN");
				break;
			case 'o':
				OutFile.open(optarg, std::ios_base::app);
				break;
			case 't':
				type.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//To have multiple output, handled by usage
	std::ostream &out = (OutFile.is_open()) ? OutFile : std::cout;

	Detector* theBox = new Detector(DetConfig, module);

	TH1D* hFlux = dynamic_cast<TH1D*> (FluxFile->Get("htotal"));
	TH1D* hXsecCC = dynamic_cast<TH1D*> (hFlux->Clone());
	TH1D* hXsecNC = dynamic_cast<TH1D*> (hFlux->Clone());
	hXsecCC->Reset("ICES");
	hXsecNC->Reset("ICES");

	double mmass = 40.0;
	double tgtXt = 1e6 * Const::Na / mmass;
	double weight = tgtXt * (theBox->WeightLAr() + 0.8 * theBox->WeightFGT());

	hFlux->Scale(1.0/pow(theBox->Zstart(), 2));

	TIter next(XsecFile->GetListOfKeys());
	TKey *kkk = dynamic_cast<TKey*> (next());
	TDirectory *dir = dynamic_cast<TDirectory*> (kkk->ReadObj());

	TGraph * gXsecCC = dynamic_cast<TGraph*> (dir->Get("tot_cc"));
	TGraph * gXsecNC = dynamic_cast<TGraph*> (dir->Get("tot_nc"));

	for (int i = 1; i < hFlux->GetNbinsX()+1; ++i)
	{
		double E = hFlux->GetBinCenter(i);
		hXsecCC->SetBinContent(i, gXsecCC->Eval(E));
		hXsecNC->SetBinContent(i, gXsecNC->Eval(E));
	}
	hXsecCC->Scale(1.0e-38);	//GENIE units for xsec is 1e-38 cm2
	hXsecNC->Scale(1.0e-38);	//GENIE units for xsec is 1e-38 cm2

	hXsecCC->Multiply(hFlux);
	hXsecNC->Multiply(hFlux);

	double tot = hXsecCC->Integral("WIDTH") + hXsecNC->Integral("WIDTH");;
	double pot = 1.0e1 * theBox->Get("Years") * theBox->Get("POT/s");

	std::cout << "CC\n";
	std::cout << "\tEvents:    " << tgtXt * 1e20 * hXsecCC->Integral("WIDTH") << "\n";
	std::cout << "\tRatio:     " << 100. * hXsecCC->Integral("WIDTH") / tot << " %\n";
	std::cout << "\tFrequency: " << 1e3 * 1e14 * weight * hXsecCC->Integral("WIDTH") << " e-3 Hz" << "\n";
	std::cout << "NC\n";
	std::cout << "\tEvents:    " << tgtXt * 1e20 * hXsecNC->Integral("WIDTH") << "\n";
	std::cout << "\tRatio:     " << 100. * hXsecNC->Integral("WIDTH") / tot << " %\n";
	std::cout << "\tFrequency: " << 1e3 * 1e14 * weight * hXsecNC->Integral("WIDTH") << " e-3 Hz" << "\n";



	out << std::setprecision(10);
	out << "LAr_" << type << "\t" << tgtXt * theBox->WeightLAr() * tot * pot << "\n";
	out << "FGT_" << type << "\t" << tgtXt * theBox->WeightFGT() * tot * pot << std::endl;

	FluxFile->Close();
	XsecFile->Close();

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
