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
	
	while((iarg = getopt_long(argc,argv, "d:l:x:f:F:o:h", longopts, &index)) != -1)
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
				OutFile.open(optarg);
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

	Detector* theBox = new Detector(DetConfig);

	TH1D* hFluxCC = dynamic_cast<TH1D*> (FluxFile->Get("htotal"));
	TH1D* hFluxNC = dynamic_cast<TH1D*> (hFluxCC->Clone());
	TH1D* hXsecCC = dynamic_cast<TH1D*> (hFluxCC->Clone());
	TH1D* hXsecNC = dynamic_cast<TH1D*> (hFluxCC->Clone());
	TH1D* hNorm = dynamic_cast<TH1D*> (hFluxCC->Clone());
	hXsecCC->Reset("ICES");
	hXsecNC->Reset("ICES");
	hNorm->Reset("ICES");

	std::cout << "at " << theBox->Zstart() << " m\n";
	std::cout << "with " << 1.0e7 * theBox->Get("Years") * theBox->Get("POT/s") << " POT\n";
	hFluxCC->Scale(1.0/pow(theBox->Zstart(), 2));
	hFluxCC->Scale(1.0e14);
	hFluxNC->Scale(1.0/pow(theBox->Zstart(), 2));
	hFluxNC->Scale(1.0e14);

	TIter next(XsecFile->GetListOfKeys());
	TKey *kkk = dynamic_cast<TKey*> (next());
	TDirectory *dir = dynamic_cast<TDirectory*> (kkk->ReadObj());

	TGraph * gXsecCC = dynamic_cast<TGraph*> (dir->Get("tot_cc"));
	TGraph * gXsecNC = dynamic_cast<TGraph*> (dir->Get("tot_nc"));

	double minX = hFluxCC->GetXaxis()->GetXmin();
	double maxX = hFluxCC->GetXaxis()->GetXmax();
	double bw = hFluxCC->GetXaxis()->GetBinWidth(1);

	for (int i = 0; i < gXsecCC->GetN(); ++i)
	{
		double x, yCC, yNC;
		gXsecCC->GetPoint(i, x, yCC);
		gXsecNC->GetPoint(i, x, yNC);
		hXsecCC->Fill(x, yCC);
		hXsecNC->Fill(x, yNC);
		hNorm->Fill(x);
	}

	hXsecCC->Divide(hNorm);
	hXsecCC->Scale(1.0e-38);	//GENIE units for xsec is 1e-38 cm2
	hXsecNC->Divide(hNorm);
	hXsecNC->Scale(1.0e-38);	//GENIE units for xsec is 1e-38 cm2

	double mmass = 40.0;
	double tgtXg = 1e6 * Const::Na / mmass;
	if (module == "FGT")
		tgtXg *= 0.8;


	hFluxCC->Multiply(hXsecCC);
	hFluxNC->Multiply(hXsecNC);

	out << std::setprecision(10);
	//std::cout << tgtXg << ", " << theBox->Weight() << ", " << hFlux->Integral("WIDTH") << std::endl;
	out << "Frequency CC: " << 1e3 * tgtXg * hFluxCC->Integral("WIDTH") << " e-3 Hz" << std::endl;
	out << "Frequency NC: " << 1e3 * tgtXg * hFluxNC->Integral("WIDTH") << " e-3 Hz" << std::endl;

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
