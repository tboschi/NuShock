#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <getopt.h>

#include "Tools.h"

#include "TFile.h"
#include "TH1.h"

void Deal(std::vector<double>& TheX, std::vector<double>& TheY, std::ifstream &InFile);
void Shift(std::vector<double>& TheX, std::vector<double>& TheY, int Ybot);
double FindFlux(std::vector<double>& TheX, std::vector<double>& TheY, double Threshold);

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"input", 	required_argument, 	0, 'i'},
		{"output", 	required_argument, 	0, 'o'},
		{"yaxis", 	required_argument, 	0, 'y'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ifstream InFile;
	TFile *OutFile;
	int Ybot = 6;


	double mmu = Const::fMMuon;
	double mpi = Const::fMKaon;
	double mel = Const::fMElectron;
	double ratio = mel*mel/(mmu*mmu) * (mpi*mpi-mel*mel)*(mpi*mpi-mel*mel)/((mpi*mpi-mmu*mmu)*(mpi*mpi-mmu*mmu));
	std::cout << "Helicity " << ratio << std::endl;
	
	while((iarg = getopt_long(argc,argv, "i:o:y:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				InFile.open(optarg);
				break;
			case 'o':
				OutFile = new TFile(optarg, "RECREATE");
				break;
			case 'y':
				Ybot = strtod(optarg, NULL);
				break;
			case 'h':
				std::cout << "Convert EPS file to ROOT" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "Eps2Root [OPTIONS]" << std::endl;
				std::cout <<"\n  -i,  --input" << std::endl;
				std::cout << "\t\tInput file is a plain file" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file is a ROOT file" << std::endl;
				std::cout <<"\n  -y,  --yaxis" << std::endl;
				std::cout << "\t\tDefine bottomline for y axix" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				return 1;
		}
	
	}

	std::string Line, Key;
	std::stringstream ssL;
	std::vector<double> All_X, All_Y, 
			    Pion_X, Pion_Y, 
			    Kaon_X, Kaon_Y, 
			    Kaon0_X, Kaon0_Y, 
			    Muon_X, Muon_Y;

	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		if (Line.find("All") != std::string::npos)
			Deal(All_X, All_Y, InFile);

		if (Line.find("Pion") != std::string::npos)
			Deal(Pion_X, Pion_Y, InFile);

		if (Line.find("Kaon") != std::string::npos)
			Deal(Kaon_X, Kaon_Y, InFile);

		if (Line.find("Ka0") != std::string::npos)
			Deal(Kaon0_X, Kaon0_Y, InFile);

		if (Line.find("Muon") != std::string::npos)
			Deal(Muon_X, Muon_Y, InFile);
	}

	Shift(All_X, All_Y, Ybot);
	Shift(Pion_X, Pion_Y, 6);
	Shift(Kaon_X, Kaon_Y, Ybot);
	Shift(Kaon0_X, Kaon0_Y, Ybot);
	Shift(Muon_X, Muon_Y, Ybot);

	TH1D* hTotal = new TH1D("htotal", "total flux", 100,0,20);
	TH1D* hPion = new TH1D("hpion", "from pion", 100,0,20);
	TH1D* hKaon = new TH1D("hkaon", "from kaon", 100,0,20);
	TH1D* hKaon0 = new TH1D("hkaon0", "from kaon0", 100,0,20);
	TH1D* hMuon = new TH1D("hmuon", "from muon", 100,0,20);
	
	double Energy;
	for (int i = 0; i < 100; ++i)
	{
		Energy = i*0.2;
		//hTotal->Fill(Energy, FindFlux(All_X, All_Y, Energy));
		hPion->Fill(Energy, FindFlux(Pion_X, Pion_Y, Energy));
		hKaon->Fill(Energy, FindFlux(Kaon_X, Kaon_Y, Energy));
		hKaon0->Fill(Energy, FindFlux(Kaon0_X, Kaon0_Y, Energy));
		hMuon->Fill(Energy, FindFlux(Muon_X, Muon_Y, Energy));
	}

	double Normalize = 1e-20 * 1.3e6*1.3e6 * 0.01*0.01; 	//nu/POT/100m/cm2/GeV  @ 1m

	//hTotal->Scale(Normalize);
	hPion->Scale(Normalize);
	hKaon->Scale(Normalize);
	hKaon0->Scale(Normalize);
	hMuon->Scale(Normalize);
/*
	double mmu = Const::fMMuon;
	double mpi = Const::fMKaon;
	double mel = Const::fMElectron;
	double ratio = mel*mel/(mmu*mmu) * (mpi*mpi-mel*mel)*(mpi*mpi-mel*mel)/((mpi*mpi-mmu*mmu)*(mpi*mpi-mmu*mmu));
	std::cout << "Helicity " << ratio << std::endl;
*/	hPion->Scale(ratio);

	hTotal->Add(hPion);
	hTotal->Add(hKaon);
	hTotal->Add(hKaon0);
	hTotal->Add(hMuon);

	OutFile->Write();

	return 0;
}
	
void Deal(std::vector<double>& TheX, std::vector<double>& TheY, std::ifstream &InFile)
{
	std::string Line;
	std::stringstream ssL;
	double PosX = 0, PosY = 0, X, Y;

	while(std::getline(InFile, Line))
	{
		if (Line[0] == '$')
			break;
		ssL.str("");
		ssL.clear();
		ssL << Line;
		ssL >> X >> Y;
		PosX += X;
		PosY += Y;
		TheX.push_back(PosX);
		TheY.push_back(PosY);
	}
}

void Shift(std::vector<double>& TheX, std::vector<double>& TheY, int Ybot)
{
	//from fit Y
	double aY = -0.499546;
	double bY = 0.00120956;
	//from fit X
	double aX = -2.49987;
	double bX = 0.00440948;
	for (int i = 0; i < TheX.size(); ++i)
	{
		TheX.at(i) = aX + bX*TheX.at(i);
		TheY.at(i) = pow(10, Ybot + aY + bY*TheY.at(i));
		if (TheY.at(i) <= pow(10, Ybot)) TheY.at(i) = 0;
	}
}

double FindFlux(std::vector<double>& TheX, std::vector<double>& TheY, double Threshold)
{
	int i = 0;
	while(TheX.at(i) <= Threshold+0.01 && i < TheX.size() )
	{
		++i;
	}
	return TheY.at(i);
}

