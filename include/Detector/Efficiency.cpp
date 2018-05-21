#include "Efficiency.h"

Efficiency::Efficiency(std::string InFile, bool Time)
{
	double Mass;
	std::string Line, SimFile, CutFile, BkgName;
	std::stringstream ssL;

	std::ifstream Input(InFile.c_str());
	while (std::getline(Input, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Mass >> SimFile >> CutFile;
		std::cout << Mass << "\t" << SimFile << "\t" << CutFile << std::endl;

		vMass.push_back(Mass);
		vSim.push_back(SimFile);
		vCut.push_back(CutFile);
	}	

	Timing = Time;
}

void Efficiency::InitFunc()	//energies from 0 to 20, 100 bin; masses from 0 to 0.5, 100 bin
{
	hhFunc = new TH2D("hhfunc", "Efficiency function", 100, 0, 20, 100, 0, 0.5);
	hCut = new TH1D("hcut", "Energy after cut", 100, 0, 20);
	hAll = new TH1D("hall", "Energy before cut", 100, 0, 20);
}

void Efficiency::LoopFile()
{
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		SetMass(vMass.at(i));
		TreeFile = new TFile(vSim.at(i).c_str(), "OPEN");
		TreeFile->cd();
	 	Data = dynamic_cast<TTree*> (TreeFile->Get("Data"));
		InitTree();
		LoadCut(vCut.at(i));

		LoopTree();
		MakeFunction();

		TreeFile->Close();
	}		
}

void Efficiency::LoopTree()
{
	hAll->Reset("ICES");
	hCut->Reset("ICES");

	for (unsigned int i = 0; i < Data->GetEntries(); ++i)
	{
		Data->GetEntry(i);

		FillAll();
		FillCut();
	}
}


void Efficiency::InitTree()
{
	Data->SetBranchAddress("Real", &Real, &b_fReal);

	Data->SetBranchAddress("E_A", &E_A, &b_fEnergyA);
	Data->SetBranchAddress("P_A", &P_A, &b_fMomentA);
	Data->SetBranchAddress("T_A", &T_A, &b_fTransvA);
	Data->SetBranchAddress("TheA", &TheA, &b_fThetaA);
	Data->SetBranchAddress("PhiA", &PhiA, &b_fPhiA);
	//Data->SetBranchAddress("M_A", &M_A, &b_fMassA);
	Data->SetBranchAddress("E_B", &E_B, &b_fEnergyB);
	Data->SetBranchAddress("P_B", &P_B, &b_fMomentB);
	Data->SetBranchAddress("T_B", &T_B, &b_fTransvB);
	Data->SetBranchAddress("TheB", &TheB, &b_fThetaB);
	Data->SetBranchAddress("PhiB", &PhiB, &b_fPhiB);
	//Data->SetBranchAddress("M_B", &M_B, &b_fMassB);
	Data->SetBranchAddress("Angle", &Angle, &b_fAngle);
	Data->SetBranchAddress("E_0", &E_0, &b_fEnergy0);
	Data->SetBranchAddress("P_0", &P_0, &b_fMoment0);
	Data->SetBranchAddress("T_0", &T_0, &b_fTransv0);
	Data->SetBranchAddress("The0", &The0, &b_fTheta0);
	//Data->SetBranchAddress("Phi0", &Phi0, &b_fPhi0);
	Data->SetBranchAddress("M_0", &M_0, &b_fMass0);
}

void Efficiency::LoadCut(std::string CutFile)
{
	Data->SetBranchStatus("*", 0);		// disable all branches
	Data->SetBranchStatus("Real", 1);	// always on
	mRef.clear();
	mCutLo.clear();
	mCutUp.clear();
	mSpecialLo.clear();
	mSpecialUp.clear();

	double Lower, Upper;
	std::string Line, BranchName;
	std::stringstream ssL;

	std::ifstream Input(CutFile.c_str());
	while (std::getline(Input, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> BranchName >> Lower >> Upper;

		if (BranchName == "E_A")
			SetMap(BranchName, &E_A, Lower, Upper);
		else if (BranchName == "P_A")
			SetMap(BranchName, &P_A, Lower, Upper);
		else if (BranchName == "T_A")
			SetMap(BranchName, &T_A, Lower, Upper);
		else if (BranchName == "TheA")
			SetMap(BranchName, &TheA, Lower, Upper);
		else if (BranchName == "PhiA")
			SetMap(BranchName, &PhiA, Lower, Upper);
		else if (BranchName == "E_B")
			SetMap(BranchName, &E_B, Lower, Upper);
		else if (BranchName == "P_B")
			SetMap(BranchName, &P_B, Lower, Upper);
		else if (BranchName == "T_B")
			SetMap(BranchName, &T_B, Lower, Upper);
		else if (BranchName == "TheB")
			SetMap(BranchName, &TheB, Lower, Upper);
		else if (BranchName == "PhiB")
			SetMap(BranchName, &PhiB, Lower, Upper);
		else if (BranchName == "E_0")
			SetMap(BranchName, &E_0, Lower, Upper);
		else if (BranchName == "P_0")
			SetMap(BranchName, &P_0, Lower, Upper);
		else if (BranchName == "T_0")
			SetMap(BranchName, &T_0, Lower, Upper);
		else if (BranchName == "The0")
			SetMap(BranchName, &The0, Lower, Upper);
		else if (BranchName == "M_0")
			SetMap(BranchName, &M_0, Lower, Upper);
		else if (BranchName == "Angle")
			SetMap(BranchName, &Angle, Lower, Upper);
		else if (BranchName == "CosAB")		//cos(PhiA-PhiB)
			SetSpecial(0, Lower, Upper);
		else if (BranchName == "aCosAB")	//abs(cos(PhiA-PhiB))
			SetSpecial(1, Lower, Upper);
		else if (BranchName == "CircAB")	//E_A*E_A+E_B*E_B
			SetSpecial(2, Lower, Upper);
		else if (BranchName == "atAB0")		//fabs(TheA-TheB)/The0
			SetSpecial(3, Lower, Upper);
		else if (BranchName == "T_A+B")		//T_A+T_B
			SetSpecial(4, Lower, Upper);
		else if (BranchName == "T_AB0")		//(T_A+T_B)/T_0
			SetSpecial(5, Lower, Upper);
		else
			std::cout << BranchName << ": branch unknown!" << std::endl;
	}
}

void Efficiency::SetMap(std::string BN, double *Address, double Lo, double Up)
{
	Data->SetBranchStatus(BN.c_str(), 1);		// activate branchname
	mRef[BN] = Address;
	mCutLo[BN] = Lo;
	mCutUp[BN] = Up;
}

void Efficiency::SetSpecial(int CutNumber, double Lo, double Up)
{
	switch (CutNumber)
	{
		case 0:		//Special cut CosAB
		case 1:		//Special cut AbsAB
			Data->SetBranchStatus("PhiA", 1);
			Data->SetBranchStatus("PhiB", 1);
			break;
		case 2:		//Special cut CircAB
			Data->SetBranchStatus("E_A", 1);
			Data->SetBranchStatus("E_B", 1);
			break;
		case 3:		//Special cut atAB0
			Data->SetBranchStatus("TheA", 1);
			Data->SetBranchStatus("TheB", 1);
			Data->SetBranchStatus("The0", 1);
			break;
		case 4:		//Special cut T_A+B
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("T_B", 1);
			break;
		case 5:		//Special cut T_AB0
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("T_B", 1);
			Data->SetBranchStatus("T_0", 1);
			break;
	}

	mSpecialLo[CutNumber] = Lo;
	mSpecialUp[CutNumber] = Up;
}

void Efficiency::FillAll()
{
	hAll->Fill(Real);
}

void Efficiency::FillCut()
{
	if (Timing && Real < 5.0)
		hCut->Fill(Real);
	else if (PassCut() && SpecialCut())
		hCut->Fill(Real);
}

bool Efficiency::PassCut()
{
	bool Ret = true;

	for (im = mRef.begin(); im != mRef.end(); ++im)
	{
		double Value = *(im->second);
		if (!(Value > mCutLo[im->first] && 
		      Value < mCutUp[im->first]))
		{
			Ret = false;
			break;
		}
	}

	return Ret;
}

bool Efficiency::SpecialCut()
{
	bool Ret = true;

	for (is = mSpecialLo.begin(); is != mSpecialLo.end(); ++is)
	{
		double Value;
		switch(is->first)
		{
			case 0:		//Special cut CosAB
				Value = cos(PhiA-PhiB);
				break;
			case 1:		//Special cut AbsAB
				Value = fabs(cos(PhiA-PhiB));
				break;
			case 2:		//Special cut CircAB
				Value = E_A*E_A + E_B*E_B;
				break;
			case 3:		//Special cut atAB0
				Value = fabs(TheA-TheB)/The0;
				break;
			case 4:		//Special cut T_A+B
				Value = T_A+T_B;
				break;
			case 5:		//Special cut T_AB0
				Value = (T_A+T_B)/T_0;
				break;
			default:
				break;
		}

		if (!(Value > is->second &&
		      Value < mSpecialUp[is->first]))
		{
			Ret = false;
			break;
		}
	}

	return Ret;
}

void Efficiency::MakeFunction()
{
	for (int Bin = 1; Bin < hAll->GetNbinsX()+1; ++Bin)
	{
		int yBin = hhFunc->GetYaxis()->FindBin(GetMass());
		double Ratio = hAll->GetBinContent(Bin) == 0 ? 1.0 : hCut->GetBinContent(Bin)/hAll->GetBinContent(Bin);

		hhFunc->SetBinContent(Bin, yBin, Ratio);
	}
}


void Efficiency::CompleteFunction()
{
	TRandom3 MT;
	TH1D *ProjX = hhFunc->ProjectionX();	//energies
	TH1D *ProjY = hhFunc->ProjectionY();	//masses

	int BinMass0, BinMass1;
	double f0, f1, M0, M1;
	double mFactor, qFactor, Mass, Eff;

	for (int Ebin = 1; Ebin < ProjX->GetNbinsX()+1; ++Ebin)
	{
		BinMass0 = ProjY->FindBin(vMass.at(0));	//0
		BinMass1 = ProjY->FindBin(vMass.at(1));	//1

		f0 = hhFunc->GetBinContent(Ebin, BinMass0);
		f1 = hhFunc->GetBinContent(Ebin, BinMass1);
		M0 = vMass.at(0);
		M1 = vMass.at(1);

		mFactor = (f1 - f0)/(M1 - M0);
		qFactor = f0 - mFactor * M0;

		for (int Mbin = 1; Mbin < BinMass1; ++Mbin)
		{
			Mass = ProjY->GetBinCenter(Mbin);
			Eff = Mass * mFactor + qFactor;
			Eff = MT.Gaus(Eff, 0.01);

			if (Eff < 0.0)
				hhFunc->SetBinContent(Ebin, Mbin, 0.0);
			else if (Eff > 1.0)
				hhFunc->SetBinContent(Ebin, Mbin, 1.0);
			else
				hhFunc->SetBinContent(Ebin, Mbin, Eff);
		}

		for (int i = 1; i < vMass.size()-1; i++)	//loop from i = 1 to i = n
		{
			BinMass0 = ProjY->FindBin(vMass.at(i));	//i
			BinMass1 = ProjY->FindBin(vMass.at(i+1));	//i+1

			f0 = hhFunc->GetBinContent(Ebin, BinMass0);
			f1 = hhFunc->GetBinContent(Ebin, BinMass1);
			M0 = vMass.at(i);
			M1 = vMass.at(i+1);

			mFactor = (f1 - f0)/(M1 - M0);
			qFactor = f0 - mFactor * M0;

			for (int Mbin = BinMass0+1; Mbin < BinMass1; ++Mbin)
			{
				Mass = ProjY->GetBinCenter(Mbin);
				Eff = Mass * mFactor + qFactor;
				Eff = MT.Gaus(Eff, 0.01);

				if (Eff < 0.0)
					hhFunc->SetBinContent(Ebin, Mbin, 0.0);
				else if (Eff > 1.0)
					hhFunc->SetBinContent(Ebin, Mbin, 1.0);
				else
					hhFunc->SetBinContent(Ebin, Mbin, Eff);
			}
		}
		//coeff known for the last bit
		for (int Mbin = BinMass1; Mbin < ProjY->GetNbinsX()+1; ++Mbin)
		{
			Mass = ProjY->GetBinCenter(Mbin);
			Eff = Mass * mFactor + qFactor;
			Eff = MT.Gaus(Eff, 0.01);

			if (Eff < 0.0)
				hhFunc->SetBinContent(Ebin, Mbin, 0.0);
			else if (Eff > 1.0)
				hhFunc->SetBinContent(Ebin, Mbin, 1.0);
			else
				hhFunc->SetBinContent(Ebin, Mbin, Eff);
		}
	}
}

TH2D *Efficiency::GetFunction()
{
	return hhFunc;
}

TH1D *Efficiency::GetAll()
{
	return hAll;
}

TH1D *Efficiency::GetCut()
{
	return hCut;
}

double Efficiency::GetMass()
{
	return dMass;
}

void Efficiency::SetMass(double X)
{
	dMass = X;
}
