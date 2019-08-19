#include "Efficiency.h"

Efficiency::Efficiency(std::string mcFile) :
	hAll(NULL),
	hCut(NULL),
	W(1)
{
	inFile = new TFile(mcFile.c_str());
	TTree *data = dynamic_cast<TTree*> (inFile->Get("Data"));

	LoadTree(data);
}

Efficiency::~Efficiency()
{
	hAll->Delete();
	hCut->Delete();
	inFile->Close();
}

void Efficiency::LoadTree(TTree *tree)
{
	if (!tree)
		return;
	Data = tree;

	Data->SetBranchAddress("E_A", &E_A, &b_fEnergyA);
	Data->SetBranchAddress("P_A", &P_A, &b_fMomentA);
	Data->SetBranchAddress("T_A", &T_A, &b_fTransvA);
	Data->SetBranchAddress("TheA", &TheA, &b_fThetaA);
	Data->SetBranchAddress("PhiA", &PhiA, &b_fPhiA);
	Data->SetBranchAddress("M_A", &M_A, &b_fMassA);
	Data->SetBranchAddress("In_A", &In_A, &b_fLengthA);
	Data->SetBranchAddress("Out_A", &Out_A, &b_fLengthoA);
	Data->SetBranchAddress("E_B", &E_B, &b_fEnergyB);
	Data->SetBranchAddress("P_B", &P_B, &b_fMomentB);
	Data->SetBranchAddress("T_B", &T_B, &b_fTransvB);
	Data->SetBranchAddress("TheB", &TheB, &b_fThetaB);
	Data->SetBranchAddress("PhiB", &PhiB, &b_fPhiB);
	Data->SetBranchAddress("M_B", &M_B, &b_fMassB);
	Data->SetBranchAddress("In_B", &In_B, &b_fLengthB);
	Data->SetBranchAddress("Out_B", &Out_B, &b_fLengthoB);
	Data->SetBranchAddress("Angle", &Angle, &b_fAngle);
	Data->SetBranchAddress("E_0", &E_0, &b_fEnergy0);
	Data->SetBranchAddress("P_0", &P_0, &b_fMoment0);
	Data->SetBranchAddress("T_0", &T_0, &b_fTransv0);
	Data->SetBranchAddress("The0", &The0, &b_fTh0ta0);
	Data->SetBranchAddress("Phi0", &Phi0, &b_fPhi0);
	Data->SetBranchAddress("M_0", &M_0, &b_fMass0);

	if (Data->GetBranch("True"))
	{
		Data->SetBranchAddress("True", &True);//, &b_fTrue);
		Data->SetBranchAddress("W",    &W);//,    &b_fW);
		Hist = &True;
	}
	else
		Hist = &E_0;
}

void Efficiency::LoadCut(std::string cutFile)
{
	Data->SetBranchStatus("*", 0);		// disable all branches
	if (Data->GetBranch("True"))
	{
		Data->SetBranchStatus("True", 1);
		Data->SetBranchStatus("W",    1);
	}
	else
		Data->SetBranchStatus("E_0", 1);

	mRef.clear();
	mCutLo.clear();
	mCutUp.clear();
	mSpecialLo.clear();
	mSpecialUp.clear();

	double lower, upper;
	std::string line, branchName;

	std::ifstream in(cutFile.c_str());
	while (std::getline(in, line))
	{
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));
		if (line.empty())
			continue;

		std::stringstream ssl(line);
		ssl >> branchName >> lower >> upper;

		if (branchName == "E_A")
			SetCut(branchName, &E_A, lower, upper);
		else if (branchName == "P_A")
			SetCut(branchName, &P_A, lower, upper);
		else if (branchName == "T_A")
			SetCut(branchName, &T_A, lower, upper);
		else if (branchName == "TheA")
			SetCut(branchName, &TheA, lower, upper);
		else if (branchName == "PhiA")
			SetCut(branchName, &PhiA, lower, upper);
		else if (branchName == "E_B")
			SetCut(branchName, &E_B, lower, upper);
		else if (branchName == "P_B")
			SetCut(branchName, &P_B, lower, upper);
		else if (branchName == "T_B")
			SetCut(branchName, &T_B, lower, upper);
		else if (branchName == "TheB")
			SetCut(branchName, &TheB, lower, upper);
		else if (branchName == "PhiB")
			SetCut(branchName, &PhiB, lower, upper);
		else if (branchName == "E_0")
			SetCut(branchName, &E_0, lower, upper);
		else if (branchName == "P_0")
			SetCut(branchName, &P_0, lower, upper);
		else if (branchName == "T_0")
			SetCut(branchName, &T_0, lower, upper);
		else if (branchName == "The0")
			SetCut(branchName, &The0, lower, upper);
		else if (branchName == "M_0")
			SetCut(branchName, &M_0, lower, upper);
		else if (branchName == "Angle")
			SetCut(branchName, &Angle, lower, upper);
		else if (branchName == "CosAB")		//cos(PhiA-PhiB)
			SetSpecial(0, lower, upper);
		else if (branchName == "aCosAB")	//abs(cos(PhiA-PhiB))
			SetSpecial(1, lower, upper);
		else if (branchName == "CircAB")	//E_A*E_A+E_B*E_B
			SetSpecial(2, lower, upper);
		else if (branchName == "atAB0")		//fabs(TheA-TheB)/The0
			SetSpecial(3, lower, upper);
		else if (branchName == "T_A+B")		//T_A+T_B
			SetSpecial(4, lower, upper);
		else if (branchName == "T_AB0")		//(T_A+T_B)/T_0
			SetSpecial(5, lower, upper);
		else if (branchName == "TTA")		//(T_A+T_B)/T_0
			SetSpecial(6, lower, upper);
		else
			std::cout << branchName << ": branch unknown!" << std::endl;
	}
}

void Efficiency::LoadSpectra(double mass)
{
	if (hAll)
		hAll->Delete();
	if (hCut)
		hCut->Delete();

	hAll = new TH1D("hall", "Energy before cut", 100, 0, 20);
	hCut = new TH1D("hcut", "Energy after cut", 100, 0, 20);

	for (int i = 0; i < Data->GetEntries(); ++i)
	{
		Data->GetEntry(i);

		hAll->Fill(*Hist, W);
		if (PassCut() && SpecialCut())
			hCut->Fill(*Hist, W);
	}

	//if (mass >= 0.0)	//store if mass is positive, for function extrapoltion
	//{
	//	if (mAll.count(mass))
	//		delete mAll[mass];
	//	if (mCut.count(mass))
	//		delete mCut[mass];

	//	mAll[mass] = hAll;
	//	mCut[mass] = hCut;
	//}
}

double Efficiency::EventsLeft()
{
	return hCut->Integral();
}

double Efficiency::ReductionFactor()
{
	return hCut->Integral() / hAll->Integral();
}

void Efficiency::SetCut(std::string BN, double *Address, double Lo, double Up)
{
	Data->SetBranchStatus(BN.c_str(), 1);		// activate branchname
	mRef[BN] = Address;
	mCutLo[BN] = Lo;
	mCutUp[BN] = Up;
}

void Efficiency::SetSpecial(int cutNumber, double Lo, double Up)
{
	switch (cutNumber)
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
		case 6:		//Special cut T_AB0
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("TheA", 1);
			break;
	}

	mSpecialLo[cutNumber] = Lo;
	mSpecialUp[cutNumber] = Up;
}

bool Efficiency::PassCut()
{
	bool ret = true;
	for (im = mRef.begin(); im != mRef.end(); ++im)
	{
		double Value = *(im->second);
		if (!(Value > mCutLo[im->first] && 
		      Value < mCutUp[im->first]))
		{
			ret = false;
			break;
		}
	}

	return ret;
}

bool Efficiency::SpecialCut()
{
	bool ret = true;
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
			case 6:		//Special cut T_AB0
				Value = T_A*TheA;
				break;
			default:
				break;
		}

		if (!(Value > is->second &&
		      Value < mSpecialUp[is->first]))
		{
			ret = false;
			break;
		}
	}

	return ret;
}

/*
TH2D* Efficiency::MakeFunction()
{
	TH2D* hhFunc = new TH2D("hhfunc", "Efficiency function", 100, 0, 20, 400, 0, 2.0);

	//loop through masses and histograms and lohad TH2
	std::map<double, TH1D*>::iterator im;
	for (im = mAll.begin(); im != mAll.end(); ++im)
	{
		for (int bin = 1; bin < im->second->GetNbinsX()+1; ++bin)
		{
			int yBin = hhFunc->GetYaxis()->FindBin(im->first);
			double frac = (im->second->GetBinContent(Bin) == 0 ? 1.0 :
					mCut[im->first]->GetBinContent(Bin) / 
					mAll[im->first]->GetBinContent(Bin);

			hhFunc->SetBinContent(bin, yBin, frac);
		}
	}

	//complete hhFunc by extrapolation

	TRandom3 MT;
	TH1D *ProjX = hhFunc->ProjectionX();	//energies
	TH1D *ProjY = hhFunc->ProjectionY();	//masses

	int BinMass0, BinMass1;
	double f0, f1, M0, M1;
	double mFactor, qFactor, Mass, Eff;

	for (int Ebin = 1; Ebin < ProjX->GetNbinsX()+1; ++Ebin)
	{
		BinMass0 = ProjY->FindFirstBin();
		BinMass1 = ProjY->FindFirstBin(0, 1, BinMass0);

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

		while (BinMass1 < ProjY->GetNbinsX())
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
*/

TH1D *Efficiency::GetAll()
{
	return hAll;
}

TH1D *Efficiency::GetCut()
{
	return hCut;
}
