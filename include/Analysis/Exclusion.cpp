#include "Analysis/Exclusion.h"

Exclusion::Exclusion(Engine* TE, Engine::Current HornType, 
		     Detector *TB,
		     std::vector<char> &vF, double Threshold) :
	TheEngine(TE),
	Horn(HornType),
	TheBox(TB),
	vFlag(vF),
	Thr(Threshold)
{
}

Exclusion::~Exclusion()
{
}

double Exclusion::Bisect(double S, double E)
{
	double fS = Function(S);
	double fE = Function(E);

	while (fabs(S - E) > 1e-3)
	{
		double M = (S + E) / 2.0;
		double fM = Function(M);

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

	return (S + E) / 2.0;
}

//the two intervals are [S : M] and [M : E]
bool Exclusion::FindInterval(double S, double &M, double E)
{
	std::list<double> ll;
	ll.push_front(S);	//this one should be less than...
	ll.push_back(E);	//...this one

	unsigned int g = 0;	//max number of splits is 50
	while ( (M = IsZero(ll)) < ll.front() && ll.size() < 100 )	//there is no zero
		Split(ll);				//keep splitting

	return (ll.size() < 100);
}

double Exclusion::IsZero(std::list<double> &ll)
{
	std::list<double>::iterator il = ll.begin();
	double f0 = Function(ll.front());
	double f1;

	for (++il; il != ll.end(); ++il)
	{
		f1 = Function(*il);

		if (f0 * f1 < 0)
			return *il;
	}

	return ll.front()-10000;	//better one needed!
}

void Exclusion::Split(std::list<double> &ll)
{
	bool Ret = false;
	std::list<double>::iterator il = ll.begin(), ip = il;
	double f0 = Function(ll.front());
	double f1;

	for (++il; il != ll.end(); ++il, ++ip)
	{
		double M = (*ip + *il) / 2.0;
		ll.insert(il, M);
	}
}

double Exclusion::Function(double lu2)
{
	SetMix(lu2);

	switch (Horn)
	{
		case Engine::FHC:
			return TheEngine->MakeSampler(TheBox, Engine::FHC) - Thr;
		case Engine::RHC:
			return TheEngine->MakeSampler(TheBox, Engine::RHC) - Thr;
		case Engine::Both:
			return TheEngine->MakeSampler(TheBox) - Thr;
		default:
			return 0.0;
	}
}

void Exclusion::SetMix(double lu2)
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
					ue = uu;
				else if (vFlag.at(f) == 'M')
					um = uu;
				else if (vFlag.at(f) == 'T')
					ut = uu;
			}
			//std::cout << ue << "\t" << um << "\t" << ut << std::endl;
			NN->SetMixings(ue, um, ut);
		}
	}
}

void Exclusion::SetThr(double T)
{
	Thr = T;
}
