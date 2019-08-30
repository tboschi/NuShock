#include "Efficiency.h"

Efficiency::Efficiency(std::string mcFile) :
	hAll(NULL),
	hCut(NULL),
	hhFunc(NULL),
	inFile(NULL),
	funcSet(false),
	W(1)
{
	MapCuts();
	LoadFile(mcFile);
}

Efficiency::~Efficiency()
{
	if (hAll)
		hAll->Delete();
	if (hCut)
		hCut->Delete();
	if (hhFunc)
		hhFunc->Delete();
	if (inFile && inFile->IsOpen())
		inFile->Close();

	hAll = NULL;
	hCut = NULL;
	hhFunc = NULL;
	inFile = NULL;
}

void Efficiency::MapCuts()
{
	mCut["E_A"]   = _E_A;
	mCut["P_A"]   = _P_A;
	mCut["T_A"]   = _T_A;
	mCut["TheA"]  = _TheA;
	mCut["E_B"]   = _E_B;
	mCut["P_B"]   = _P_B;
	mCut["T_B"]   = _T_B;
	mCut["TheB"]  = _TheB;
	mCut["Angle"] = _Angle;
	mCut["E_0"]   = _E_0;
	mCut["P_0"]   = _P_0;
	mCut["T_0"]   = _T_0;
	mCut["The0"]  = _The0;
	mCut["M_0"]   = _M_0;
	                     
	mCut["CosAB"]  = _CosAB;
	mCut["aCosAB"] = _aCosAB;
	mCut["CircAB"] = _CircAB;
	mCut["atAB0"]  = _atAB0;
	mCut["T_A+B"]  = _T_AB;
	mCut["T_AB0"]  = _T_AB0;
	mCut["TTA"]    = _TTA;
	mCut["EAB"]    = _EAB;
	mCut["E0Ang"]  = _E0Ang;
}

void Efficiency::LoadFile(std::string mcFile)
{
	if (!mcFile.empty())
	{
		if (inFile && inFile->IsOpen())
			inFile->Close();

		inFile = new TFile(mcFile.c_str());
		TTree *data = dynamic_cast<TTree*> (inFile->Get("Data"));
		LoadTree(data);
	}
}

std::vector<std::string> Efficiency::AvailableCuts(std::string name)
{
	std::vector<std::string> vCut;
	if (name == "nPI0")
	{
		vCut.push_back("E_0");
		vCut.push_back("T_0");
		vCut.push_back("The0");
		vCut.push_back("E_A");
		vCut.push_back("TheA");
	}
	else if (name == "nEE" || name == "nMM")
	{
		vCut.push_back("E_A");
		vCut.push_back("E_0");
		vCut.push_back("M_0");
		vCut.push_back("T_A");
		vCut.push_back("T_0");
		vCut.push_back("TheA");
		vCut.push_back("The0");
		vCut.push_back("T_AB0");
		vCut.push_back("Angle");
	}
	else if (name == "nEM")
	{
		vCut.push_back("E_A");
		vCut.push_back("E_B");	//0.5
		vCut.push_back("E_0");
		vCut.push_back("M_0");
		vCut.push_back("T_A");
		vCut.push_back("T_B");
		vCut.push_back("T_0");
		vCut.push_back("TheA");
		vCut.push_back("TheB");
		vCut.push_back("The0");
		vCut.push_back("T_AB0");
		vCut.push_back("Angle");
	}
	else if (name == "EPI" || name == "MPI")
	{
		vCut.push_back("E_A");
		vCut.push_back("E_B");	//0.5
		vCut.push_back("E_0");
		vCut.push_back("M_0");
		vCut.push_back("T_A");
		vCut.push_back("T_B");
		vCut.push_back("T_0");
		vCut.push_back("TheA");
		vCut.push_back("TheB");
		vCut.push_back("The0");
		vCut.push_back("Angle");
		vCut.push_back("CosAB");
		vCut.push_back("T_AB0");
	}

	return vCut;

}

void Efficiency::LoadTree(TTree *tree)
{
	if (!tree)
		return;
	Data = tree;

	Data->SetBranchAddress("PdgA", &PdgA, &b_fPdgA);
	Data->SetBranchAddress("E_A", &E_A, &b_fEnergyA);
	Data->SetBranchAddress("P_A", &P_A, &b_fMomentA);
	Data->SetBranchAddress("T_A", &T_A, &b_fTransvA);
	Data->SetBranchAddress("TheA", &TheA, &b_fThetaA);
	Data->SetBranchAddress("PhiA", &PhiA, &b_fPhiA);
	Data->SetBranchAddress("PdgB", &PdgB, &b_fPdgB);
	Data->SetBranchAddress("E_B", &E_B, &b_fEnergyB);
	Data->SetBranchAddress("P_B", &P_B, &b_fMomentB);
	Data->SetBranchAddress("T_B", &T_B, &b_fTransvB);
	Data->SetBranchAddress("TheB", &TheB, &b_fThetaB);
	Data->SetBranchAddress("PhiB", &PhiB, &b_fPhiB);
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

//CL regulates the threshold, it varies between 0 and 1
//until background increases
void Efficiency::GetCutLimits(std::string name, double &cLo, double &cUp)
{
	if (mRef.count(name))
	{
		cLo = mLower[name];
		cUp = mUpper[name];
	}
}

void Efficiency::FindCut(std::string name, double &cLo, double &cUp, double CL, double mass)
{
	SetCut(name);

	TH1D* hist = LoadCutSpectrum(name);
	double eff = hist->GetEntries() * 1e-4;

	int sL = 1, sR = hist->GetNbinsX();
	if (mass > 0 && name == "M_0")
		FindRange(hist, sL, sR, hist->FindBin(mass), 0.8*(1-CL)+CL);
	else if (name == "E_A" || name == "E_B")
		FindRange(hist, sL, sR, hist->GetNbinsX(), 0.70*(1-CL)+CL);
	else if (name == "E_0")	// || name == "Angle")
		FindRange(hist, sL, sR, hist->GetNbinsX(), 0.9*(1-CL)+CL);
	else if (name == "T_AB0")
		FindRange(hist, sL, sR, 1, 0.90*(1-CL)+CL);
	else if (name == "T_A" || name == "T_B" || name == "T_0")
	{
		int s0 = hist->GetMaximumBin();
		double thr = hist->GetBinContent(s0) * 0.05 * (1-CL);	//5% of max
		sL = hist->FindFirstBinAbove(thr);
		sR = hist->FindLastBinAbove(thr);
	}
	else if (name == "TheA" || name == "TheB" ||
		 name == "The0" || name == "Angle" || name == "CosAB")
	{
		int s0 = hist->GetMaximumBin();
		double thr = hist->GetBinContent(s0) * 0.10 * (1-CL);	//10.% of max
		sL = hist->FindFirstBinAbove(thr);
		sR = hist->FindLastBinAbove(thr);
	}
	//else
		//FindRange(hist, sL, sR, hist->GetMaximumBin(), 0.99);
	

	cLo = hist->GetBinLowEdge(sL);
	cUp = hist->GetBinLowEdge(sR+1);

	mLower[name] = cLo;
	mUpper[name] = cUp;

	delete hist;
}

void Efficiency::FindRange(TH1D *hist, int &sL, int &sR, int s0, double CL)
{
	sL = s0-1, sR = s0+1;
	double tot = hist->Integral();		// == num of entries
	double sum = hist->GetBinContent(s0) / tot;

	while ( sum < CL &&
	       ( sL > 0 || sR < hist->GetNbinsX()) ) //reached both ends
	{
		if (sL > 0)
		{
			sum += hist->GetBinContent(sL) / tot;
			--sL;
		}

		if (sR < hist->GetNbinsX())
		{
			sum += hist->GetBinContent(sR) / tot;
			++sR;
		}
	}

	if (sL < 1)
		sL = 1;
	if (sR > hist->GetNbinsX())
		sR = hist->GetNbinsX();
}

int Efficiency::GetEntries()
{
	return Data->GetEntries();
}

void Efficiency::GetEntry(int i)
{
	Data->GetEntry(i);

	CosAB  = cos(PhiA-PhiB);
	aCosAB = fabs(cos(PhiA-PhiB));
	CircAB = E_A*E_A + E_B*E_B;
	atAB0  = fabs(TheA-TheB)/The0;
	T_AB   = T_A+T_B;
	T_AB0  = T_0 - std::abs(T_A-T_B);
	TTA    = T_A*TheA;
	EAB    = E_A / E_B;
	E0Ang  = E_0 * Angle;
}

TH1D* Efficiency::LoadCutSpectrum(std::string name)
{
	GetEntry(0);
	double minVal = *mRef[name], maxVal = *mRef[name];

	for (int i = 1; i < Data->GetEntries(); ++i)
	{
		GetEntry(i);
		if (!PassCut(name))
			continue;

		if (minVal > *mRef[name])
			minVal = *mRef[name];
		if (maxVal < *mRef[name])
			maxVal = *mRef[name];
	}

	std::string hname = "hist_" + name;
	TH1D* hvar = new TH1D(hname.c_str(), name.c_str(), 100, minVal, maxVal);

	for (int i = 0; i < Data->GetEntries(); ++i)
	{
		GetEntry(i);
		if (PassCut(name))
			hvar->Fill(*mRef[name]);
	}

	return hvar;
}

void Efficiency::LoadCut(std::string cutFile)
{
	Data->SetBranchStatus("*", 0);		// disable all branches
	Data->SetBranchStatus("PdgA", 1);
	Data->SetBranchStatus("PdgB", 1);
	if (Data->GetBranch("True"))
	{
		Data->SetBranchStatus("True", 1);
		Data->SetBranchStatus("W",    1);
	}
	else
		Data->SetBranchStatus("E_0", 1);

	Reset();

	CardDealer cf(cutFile);
	std::map<std::string, std::vector<double> > mName;
	std::map<std::string, std::vector<double> >::iterator im;
	if (!cf.GetAll(mName))
	{
		std::cerr << "invalid file " << cutFile << std::endl;
		return;
	}

	for (im = mName.begin(); im != mName.end(); ++im)
		if (im->second.size() >= 2)
			SetCut(im->first, im->second.at(0), im->second.at(1));

	/*
	double lower, upper;
	std::string line, name;

	std::ifstream in(cutFile.c_str());
	while (std::getline(in, line))
	{
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));
		if (line.empty())
			continue;

		std::stringstream ssl(line);
		ssl >> name >> lower >> upper;

		SetCut(name, lower, upper);
	}
	*/
}

void Efficiency::Reset()
{
	mRef.clear();
	mLower.clear();
	mUpper.clear();
	if (hAll)
	{
		hAll->Delete();
		hAll = NULL;
	}
	if (hCut)
	{
		hCut->Delete();
		hCut = NULL;
	}
}

void Efficiency::SetCut(std::string name, double lower, double upper)
{
	Data->SetBranchStatus("PdgA", 1);
	Data->SetBranchStatus("PdgB", 1);

	if (Data->GetBranch(name.c_str()))
		Data->SetBranchStatus(name.c_str(), 1);

	double *var;
	switch (mCut[name])
	{
		case Cut::_E_A:
			var = &E_A;
			break;
		case Cut::_P_A:
			var = &P_A;
			break;
		case Cut::_T_A:
			var = &T_A;
			break;
		case Cut::_TheA:
			var = &TheA;
			break;
		case Cut::_E_B:
			var = &E_B;
			break;
		case Cut::_P_B:
			var = &P_B;
			break;
		case Cut::_T_B:
			var = &T_B;
			break;
		case Cut::_TheB:
			var = &TheB;
			break;
		case Cut::_Angle:
			var = &Angle;
			break;
		case Cut::_E_0:
			var = &E_0;
			break;
		case Cut::_P_0:
			var = &P_0;
			break;
		case Cut::_T_0:
			var = &T_0;
			break;
		case Cut::_The0:
			var = &The0;
			break;
		case Cut::_M_0:
			var = &M_0;
			break;
		case Cut::_CosAB:
			Data->SetBranchStatus("PhiA", 1);
			Data->SetBranchStatus("PhiB", 1);
			var = &CosAB;
			break;
		case Cut::_aCosAB:
			Data->SetBranchStatus("PhiA", 1);
			Data->SetBranchStatus("PhiB", 1);
			var = &aCosAB;
			break;
		case Cut::_CircAB:
			Data->SetBranchStatus("E_A", 1);
			Data->SetBranchStatus("E_B", 1);
			var = &CircAB;
			break;
		case Cut::_atAB0:
			Data->SetBranchStatus("TheA", 1);
			Data->SetBranchStatus("TheB", 1);
			Data->SetBranchStatus("The0", 1);
			var = &atAB0;
			break;
		case Cut::_T_AB:
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("T_B", 1);
			var = &T_AB;
			break;
		case Cut::_T_AB0:
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("T_B", 1);
			Data->SetBranchStatus("T_0", 1);
			var = &T_AB0;
			break;
		case Cut::_TTA:
			Data->SetBranchStatus("T_A", 1);
			Data->SetBranchStatus("TheA", 1);
			var = &TTA;
			break;
		case Cut::_E0Ang:
			Data->SetBranchStatus("E_0", 1);
			Data->SetBranchStatus("Angle", 1);
			var = &E0Ang;
			break;
		default:
			std::cout << "Cut unknown: " << name << std::endl;
	}

	mRef[name] = var;
	mLower[name] = lower;
	mUpper[name] = upper;
}

void Efficiency::ApplyCut(double mass)
{
	if (hAll)
		hAll->Delete();
	if (hCut)
		hCut->Delete();

	hAll = new TH1D("hall", "Energy before cut", 100, 0, 20);
	hCut = new TH1D("hcut", "Energy after cut", 100, 0, 20);
	hAll->SetDirectory(0);
	hCut->SetDirectory(0);

	for (int i = 0; i < Data->GetEntries(); ++i)
	{
		GetEntry(i);

		hAll->Fill(*Hist, W);
		if (PassCut())
			hCut->Fill(*Hist, W);
	}

	if (funcSet)
	{
		int y = hhFunc->GetYaxis()->FindBin(mass);
		for (int x = 1; x < hAll->GetNbinsX()+1; ++x)
		{
			double frac = hAll->GetBinContent(x) == 0 ? 1.0 :
				      hCut->GetBinContent(x) / hAll->GetBinContent(x);

			hhFunc->SetBinContent(x, y, frac);
		}
	}
}

int Efficiency::EntriesLeft()
{
	return hCut->GetEntries();
}

double Efficiency::EventsLeft()
{
	return hCut->Integral();
}

double Efficiency::ReductionFactor()
{
	TH1D* hRatio = dynamic_cast<TH1D*>(hCut->Clone("reduction"));
	hRatio->Reset("ICES");
	for (int b = 1; b < hRatio->GetNbinsX()+1; ++b)
	{
		double frac = hAll->GetBinContent(b) == 0 ? 1.0 :
			hCut->GetBinContent(b) / hAll->GetBinContent(b);

		hRatio->SetBinContent(b, frac);
	}

	return hRatio->Integral() / hRatio->GetNbinsX();
}

bool Efficiency::PassCut(std::string name)
{
	for (im = mRef.begin(); im != mRef.end(); ++im)
	{
		if (im->first == name)
			continue;

		double val = *(im->second);
		if (!(val >= mLower[im->first] && 
		      val <= mUpper[im->first]))
			return false;
	}

	return true;
}

TH1D *Efficiency::GetAll()
{
	return hAll;
}

TH1D *Efficiency::GetCut()
{
	return hCut;
}

void Efficiency::MakeFunction()
{
	if (hhFunc)
		hhFunc->Delete();

	hhFunc = new TH2D("hhfunc", "Efficiency function", 100, 0, 20, 400, 0, 2.0);
	hhFunc->SetDirectory(0);
	funcSet = true;
}

int Efficiency::FindFirstBin(TH1D* hist, double thr, int start, int end)
{
	if (start < 0)
		start = 0;
	if (end < 0 || end > hist->GetNbinsX())
		end = hist->GetNbinsX();

	for (int b = start+1; b <= end; ++b)
		if (hist->GetBinContent(b) > thr)
			return b;

	return -1;
}

TH2D* Efficiency::CompleteFunction()
{
	//hhFunc must be filled
	if (!hhFunc->GetEntries())
		return 0;

	//complete hhFunc by extrapolation
	TRandom3 MT;
	TH1D *projX = hhFunc->ProjectionX();	//energies
	TH1D *projY = hhFunc->ProjectionY();	//masses

	double f0, f1, m0, m1;
	double mF, qF, mass, eff;

	for (int eBin = 1; eBin < projX->GetNbinsX()+1; ++eBin)
	{
		int mBin0 = FindFirstBin(projY);
		int mBin1 = FindFirstBin(projY, 0.0, mBin0);
		int mStart = 1, mEnd = mBin1;

		while (mBin0 > 0 && mBin0 < projY->GetNbinsX()+1 &&
		       mBin1 > 0 && mBin1 < projY->GetNbinsX()+1)
		{
			double f0 = hhFunc->GetBinContent(eBin, mBin0);
			double f1 = hhFunc->GetBinContent(eBin, mBin1);
			double m0 = projY->GetBinCenter(mBin0);
			double m1 = projY->GetBinCenter(mBin1);

			double mF = (f1 - f0)/(m1 - m0);
			double qF = f0 - mF * m0;

			for (int mBin = mStart; mBin < mEnd; ++mBin)
			{
				mass = projY->GetBinLowEdge(mBin);
				double eff = mass * mF + qF;
				eff = MT.Gaus(eff, 0.01);

				if (eff < 0.0)
					hhFunc->SetBinContent(eBin, mBin, 0.0);
				else if (eff > 1.0)
					hhFunc->SetBinContent(eBin, mBin, 1.0);
				else
					hhFunc->SetBinContent(eBin, mBin, eff);
			}

			mBin0 = mBin1;
			mBin1 = FindFirstBin(projY, 0, mBin0);

			mStart = mBin0;
			if (FindFirstBin(projY, 0, mBin1) < 0)
				mEnd = projY->GetNbinsX()+1;
			else
				mEnd = mBin1;
		}
	}

	funcSet = false;
	return hhFunc;
}
