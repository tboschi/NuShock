#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Tools.h"

double Exponent(double Tau, double Energy, double Mass);
int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"root", 	required_argument, 	0, 'r'},
		{"output", 	required_argument, 	0, 'o'},
		{"momentum", 	required_argument, 	0, 'p'},
		{"energy", 	required_argument, 	0, 'e'},
		{"lambda", 	required_argument, 	0, 'l'},
		{"length", 	required_argument, 	0, 'L'},
		{"radius", 	required_argument, 	0, 'R'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	
//Initialize variables
	
	double PionP, NeutE, Lambda, Length, Radius;
	TFile *OutF;

	while((iarg = getopt_long(argc,argv, "o:r:p:e:l:L:R:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'o':
				OutFile.open(optarg);
				break;
			case 'r':
				OutF = new TFile(optarg, "RECREATE");
				break;
			case 'p':
				PionP = strtod(optarg, NULL);
				break;
			case 'e':
				NeutE = strtod(optarg, NULL);
				break;
			case 'l':
				Lambda = strtod(optarg, NULL);
				break;
			case 'L':
				Length = strtod(optarg, NULL);
				break;
			case 'R':
				Radius = strtod(optarg, NULL);
				break;
			case 'h':
				std::cout << "Description" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "name [OPTIONS]" << std::endl;
				std::cout <<"\n  -r,  --root" << std::endl;
				std::cout << "\t\tROOT file" << std::endl;
				std::cout <<"\n  -p,  --momentum" << std::endl;
				std::cout << "\t\tPion momentum" << std::endl;
				std::cout <<"\n  -e,  --energy" << std::endl;
				std::cout << "\t\tNeutrino energy" << std::endl;
				std::cout <<"\n  -l,  --lambda" << std::endl;
				std::cout << "\t\tDecay pipe length" << std::endl;
				std::cout <<"\n  -L,  --length" << std::endl;
				std::cout << "\t\tDetector distance" << std::endl;
				return 1;
			default:
				return 1;
		}
	}

	TRandom3 *RanGen = new TRandom3(0);
	/*
	TH1D * PionSpectrum = new TH1D("pion", "pion spectrum", 500, 0, 5);
	TH1D * pPion = new TH1D("ppion", "pion spectrum", 500, 0, 5);
	TH1D * KaonSpectrum = new TH1D("kaon", "kaon spectrum", 500, 0, 5);
	TH1D * pKaon = new TH1D("pkaon", "kaon spectrum", 500, 0, 5);
	TH1D * piNeutSpectrum = new TH1D("pineut", "pionneutrino spectrum", 500, 0, 5);
	TH1D * kaNeutSpectrum = new TH1D("kaneut", "kaonneutrino spectrum", 500, 0, 5);
	TH1D * SterileSpectrum = new TH1D("sterile", "neutrino spectrum", 500, 0, 5);
	*/

	double EnergyA;
	double MomentA;
	double ThetaA;
	double PhiA;
	double px; 
	double py;  
	double pz; 

	double EnergyB;
	double MomentB;
	double ThetaB;
	double PhiB;

	double Dist;
	double ToDetect;
	double Center;
	bool InDetect;

	double Mp = Const::fMPion; 	//Pion mass
	double Mk = Const::fMKaon; 	//Kaon mass
	double Me = Const::fMElectron; 	//Elec mass
	double Mm = Const::fMMuon; 	//Elec mass
	double Mn = 0.0;		//Neutrino mass

	TTree *tNeutrino = new TTree("tNeutrino", "Neutrinos");
	tNeutrino->Branch("ThetaA", &ThetaA, "fThetaA/D");
	tNeutrino->Branch("PhiA", &PhiA, "fPhiA/D");
	tNeutrino->Branch("EnergyA", &EnergyA, "fEnergyA/D");
	tNeutrino->Branch("MomentA", &MomentA, "fMomentA/D");
	tNeutrino->Branch("ThetaB", &ThetaB, "fThetaB/D");
	tNeutrino->Branch("PhiB", &PhiB, "fPhiB/D");
	tNeutrino->Branch("EnergyB", &EnergyB, "fEnergyB/D");
	tNeutrino->Branch("MomentB", &MomentB, "fMomentB/D");
	tNeutrino->Branch("Dist", &Dist, "fDist/D");
	tNeutrino->Branch("ToDetect", &ToDetect, "fToDetect/D");
	tNeutrino->Branch("Center", &Center, "fCenter/D");
	tNeutrino->Branch("InDetect", &InDetect, "bInDetect/O");

	int Reps = 10000;
	std::vector<std::vector<double> > vMatrix;
	for (PionP = 0.0; PionP < 10; PionP += 0.1)
	{	
		TLorentzVector Pion_vec(0, 0, PionP, sqrt(PionP*PionP + Mp*Mp));
		TVector3 bbb(Pion_vec.BoostVector());

		std::vector<double> vRow;
		for (double NeutE = 0.0; NeutE < 10; NeutE += 0.1)
		{
			EnergyA = (Mp*Mp + Mn*Mn - Mm*Mm)/(2*Mp);	//pion into mu and n
			MomentA = sqrt(EnergyA*EnergyA - Mn*Mn);

			int Count = 0;
			for (int i = 0; i < Reps; ++i)
			{
				//rest frame
				ThetaA = RanGen->Uniform(-Const::fPi, Const::fPi);
				PhiA = RanGen->Uniform(0, 2 * Const::fPi);
				px = MomentA * cos(PhiA) * sin(ThetaA);
				py = MomentA * sin(PhiA) * sin(ThetaA);
				pz = MomentA * cos(ThetaA);
				TLorentzVector N_A(px, py, pz, EnergyA);
			
				//lab frame
				TLorentzVector N_B(N_A);
				N_B.Boost(bbb);
				EnergyB = N_B.E();
				MomentB = N_B.P();
				ThetaB = N_B.Theta();
				PhiB = N_B.Phi();
		
				double Exp = log(1-RanGen->Rndm())/Exponent(26.03e-9, Pion_vec.E(), Pion_vec.M());
				Dist = Exp / Const::fM2GeV;
				ToDetect = Length-Dist;
				Center = ToDetect*tan(ThetaB);
		
				if (ThetaB < Const::fPi/2 && Dist < Lambda && Center < Radius)
				{
					InDetect = true;
					if (abs(EnergyB - NeutE) < 0.05)
						++Count;
				}
				else InDetect = false;
		
				//tNeutrino->Fill();
			}
			vRow.push_back(double(Count) / Reps);
		}
		vMatrix.push_back(vRow);
	}

	for (int i = 0; i < vMatrix.size(); ++i)
	{
		for (int j = 0; j < vMatrix.at(i).size(); ++j)
		{
			OutFile << i*0.1 << "\t" << j*0.1 << "\t";
			OutFile << vMatrix.at(i).at(j) << std::endl;
		}
	}

	//PionSpectrum->Write();

	//tNeutrino->Write();

	//OutF->Close();

	return Dist;
}
	
double Exponent(double Tau, double Energy, double Mass)
{
	double Gamma = Const::fhBar/Tau;
	double Lorentz = Mass/sqrt(Energy*Energy - Mass*Mass);
	return -Lorentz*Gamma;
}	
