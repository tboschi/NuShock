#include "Analysis/ExclusionLNV.h"

ExclusionLNV::ExclusionLNV(Engine* TE, Detector *TB, std::vector<char> &vF) :
	TheEngine(TE),
	TheBox(TB),
	vFlag(vF)
{
	///////edit for LNV
	std::ifstream inFile("data/FeldmanCousinLUT.dat");
	std::string line;
	std::stringstream parse;
	int bb;
	double ss;
	while (std::getline(inFile, line))
	{
		parse.clear();
		parse.str("");

		parse << line;
		parse >> bb >> ss;

		//vbb.push_back(bb);
		vss.push_back(ss);
	}	
	inFile.close();
	///////loaded LUT for F&C background
}

ExclusionLNV::~ExclusionLNV()
{
}

	///////edit for LNV
double ExclusionLNV::fcLUT(double bb)
{
	bb = std::abs(bb);
	if (bb < vss.size())
		return vss.at(int(bb));
	else
		return vss.back();
}
	///////LUT function

double ExclusionLNV::Bisect(double S, double E, double &nFunction)
{
	double fS = ZeroCross(S);
	double fE = ZeroCross(E);

	while (fabs(S - E) > 1e-3)
	{
		double M = (S + E) / 2.0;
		double fM = ZeroCross(M);

		if (fM == 0)
		{
			S = M;
			E = M;
		}
		else if (fS * fM < 0)
		{
			E = M;
			fE = fM;
		}
		else if (fE * fM < 0)
		{
			S = M;
			fS = fM;
		}
	}

	//std::cout << "back " << ithr << "\tevts " << fS << "\t" << fE << std::endl;
	nFunction = Function((S+E)/2.0);
	return (S + E) / 2.0;
}

//the two intervals are [S : M] and [M : E]
bool ExclusionLNV::FindInterval(double S, double &M, double E)
{
	std::list<double> ll;
	ll.push_front(S);	//this one should be less than...
	ll.push_back(E);	//...this one

	unsigned int g = 0;	//max number of splits is 50
	while ( (M = IsZero(ll)) < ll.front() && ll.size() < 100 )	//there is no zero
		Split(ll);				//keep splitting

	return (ll.size() < 100);
}

double ExclusionLNV::IsZero(std::list<double> &ll)
{
	std::list<double>::iterator il = ll.begin();
	double f0 = ZeroCross(ll.front());
	double f1;

	for (++il; il != ll.end(); ++il)
	{
		f1 = ZeroCross(*il);

		if (f0 * f1 < 0)
			return *il;
	}

	return ll.front()-10000;	//better one needed!
}

void ExclusionLNV::Split(std::list<double> &ll)
{
	bool Ret = false;
	std::list<double>::iterator il = ll.begin(), ip = il;
	//double f0 = ZeroCross(ll.front());
	//double f1;

	for (++il; il != ll.end(); ++il, ++ip)
	{
		double M = (*ip + *il) / 2.0;
		ip = ll.insert(il, M);
	}
}

double ExclusionLNV::Function(double lu2)
{
	SetMix(lu2);

	///////edit for LNV
	//ithr = fcThreshold;
	double numEvents = (TheEngine->MakeSampler(TheBox, Engine::FHC, 0) + 
			    TheEngine->MakeSampler(TheBox, Engine::RHC, 0) ) / 2.0;
	//std::cout << "back " << intrinsicBack << "\t" << numEvents << "\tthr " << fcThreshold << std::endl;
	////then use fcThreshold instead of thr

	return numEvents;
}

double ExclusionLNV::ZeroCross(double lu2)
{
	double fcThreshold = fcLUT(TheEngine->MakeSampler(TheBox, Engine::RHC, 1));
	return Function(lu2) - fcThreshold;
}

void ExclusionLNV::SetMix(double lu2)
{
	double uu = pow(10, 0.5 * lu2);

	for (unsigned int h = 0; h < 2; ++h)
	{
		Neutrino *NN;
		Engine::Current Horn;
		switch(h)
		{
			case 0:
				Horn = Engine::FHC;
				break;
			case 1:
				Horn = Engine::RHC;
				break;
		}

		for (unsigned int n = 0; n < TheEngine->vNeutrino(Horn); ++n)
		{
			NN = TheEngine->vNeutrino(Horn, n);
			double ue = 0.0, um = 0.0, ut = 0.0;
			for (unsigned int f = 0; f < vFlag.size(); ++f)
			{
				if (vFlag.at(f) == 'E')
				{
					ue = uu;
					um = uu/10.;
					ut = uu/10.;
				}
				else if (vFlag.at(f) == 'M')
				{
					ue = uu/10.;
					um = uu;
					ut = uu/10.;
				}
				else if (vFlag.at(f) == 'T')
				{
					ue = uu/10.;
					um = uu/10.;
					ut = uu;
				}
			}
			NN->SetMixings(ue, um, ut);
		}
	}
}
