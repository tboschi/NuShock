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
	std::string DetConfig;
	TFile *XsecFile, *FluxFile;
	std::ofstream OutFile;
	std::string type;
	
	while((iarg = getopt_long(argc,argv, "d:x:f:t:o:h", longopts, &index)) != -1)
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

	Detector * TheBox = new Detector(DetConfig);

	TH1D* hFlux = dynamic_cast<TH1D*> (FluxFile->Get("htotal"));
	hFlux->Scale(1.0/pow(TheBox->Get("Baseline"), 2));
	//hFlux->Scale(TheBox->Get("POT/s"));
	//hFlux->Scale(TheBox->Area()*1.0e4);

	TIter Next(XsecFile->GetListOfKeys());
	TKey *Kkk = dynamic_cast<TKey*> (Next());
	TDirectory *Dir = dynamic_cast<TDirectory*> (Kkk->ReadObj());
	TGraph * gXsecCC = dynamic_cast<TGraph*> (Dir->Get("tot_cc"));
	TGraph * gXsecNC = dynamic_cast<TGraph*> (Dir->Get("tot_nc"));

	//TFile* OutTemp = new TFile("OutTemp.root", "RECREATE");
	TH1D * hXsecCC = new TH1D("hxseccc", "hxseccc", hFlux->GetNbinsX(), 0, 20);
	TH1D * hXsecNC = new TH1D("hxsecnc", "hxsecnc", hFlux->GetNbinsX(), 0, 20);

	double x = 0, y = 0;
	unsigned int j = 0;
	for (unsigned int i = 1; i < hFlux->GetNbinsX()+1; ++i)
	{
		double sumCC = 0, sumNC = 0;
		unsigned int NCC = 0, NNC = 0;
		while (x < hFlux->GetBinCenter(i) + 0.5*20.0/hFlux->GetNbinsX())
		{
			gXsecCC->GetPoint(j, x, y);
			sumCC += y;
			++NCC;
			gXsecNC->GetPoint(j, x, y);
			sumNC += y;
			++NNC;
			++j;
		}
		hXsecCC->SetBinContent(i, 1.0e-38*sumCC/NCC);
		hXsecNC->SetBinContent(i, 1.0e-38*sumNC/NNC);
	}
	//hXsecCC->Write();
	//hXsecNC->Write();

	double MolarMass = 40.0;
	/*
	switch (TheBox->GetMaterial("InTarget"))
	{
		case Detector::LAr:
		case Detector::GasAr:
			MolarMass = 40.0;	//grams
			break;
		case Detector::Fe:
			MolarMass = 56.0;	//grams
			break;
	}
	*/

	//double NtargetTotLAr = Const::Na * TheBox->Get("WeightLAr") * 1e6 / MolarMass;	//number of targets total tons
	//double NtargetTotFGT = Const::Na * TheBox->Get("WeightFGT") * 1e6 / MolarMass;	//number of targets total tons
	double NtargetPer = Const::Na * 1e6 / MolarMass;				//number of targets per ton

	double RateCC = 0, RateNC = 0;
	double Npot = 0;
	for (unsigned int i = 1; i < hFlux->GetNbinsX()+1; ++i)
	{
		Npot += hFlux->GetBinContent(i) * 0.2;
		RateCC += hFlux->GetBinContent(i) * hXsecCC->GetBinContent(i) * 0.2;
		RateNC += hFlux->GetBinContent(i) * hXsecNC->GetBinContent(i) * 0.2;
		std::cout << i << " : " << hFlux->GetBinContent(i) << "\t" << hXsecCC->GetBinContent(i) << std::endl;
	}
	std::cout << "Rate CC " << RateCC << "\tRate NC " << RateNC << std::endl;
	//RateCC *= TheBox->GetElement("Fiducial");		//fiducial volume
	//RateNC *= TheBox->GetElement("Fiducial");		//fiducial volume

	FluxFile->Close();
	XsecFile->Close();

	double Ysec = 1.0e7 * TheBox->Get("Years");
	double pot = 1.0e7 * TheBox->Get("Years") * TheBox->Get("POT/s");
	double weight = TheBox->Get("WeightLAr") + TheBox->Get("WeightFGT");
	std::cout << "Scale " << NtargetPer*TheBox->Get("WeightLAr")*(RateCC+RateNC)*pot/1e6 << std::endl;
	std::cout << "Scale " << NtargetPer*TheBox->Get("WeightFGT")*(RateCC+RateNC)*pot/1e6 << std::endl;
	out << std::setprecision(10);
	out << "LAr_" << type << "\t" << NtargetPer*TheBox->Get("WeightLAr")*(RateCC+RateNC)*pot/1e6 << std::endl;
	out << "FGT_" << type << "\t" << NtargetPer*TheBox->Get("WeightFGT")*(RateCC+RateNC)*pot/1e6 << std::endl;
	//Out << "Tot numb of neutrinos per 1e20 POT is " << Npot*1e20 << std::endl;
	//Out << "Total number of CC events per 1e20 POT per ton is " << NtargetPer*RateCC*1e20 << std::endl;
	//Out << "Total number of NC events per 1e20 POT per ton is " << NtargetPer*RateNC*1e20 << std::endl;
	//Out << "Percentage of CC over total " << RateCC/(RateNC+RateCC) << std::endl;
	//Out << "Frequency of CC events is " << 1000.*NtargetPer*RateCC*TheBox->Get("POT/s") << " e-3 Hz on average" << std::endl;
	//Out << "Frequency of NC events is " << 1000.*NtargetPer*RateNC*TheBox->Get("POT/s") << " e-3 Hz on average" << std::endl;
	//Out << "Percentage of NC over total " << RateNC/(RateNC+RateCC) << std::endl;
	//Out << "Frequency of total events is " << NtargetPer*(RateNC+RateCC)*TheBox->Get("POT/s")*Ysec << " Hz on average" << std::endl;

	//OutTemp->Close();
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
