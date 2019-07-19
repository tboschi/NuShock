#include "Exclusion.h"

Exclusion::Exclusion(Engine* TE, Detector *TB, Engine::Current type, 
		     bool ue, bool um, bool ut, double thr, double mod) :
	engine(TE),
	box(TB),
	horn(type),
	UeFlag(ue),
	UmFlag(um),
	UtFlag(ut),
	threshold(thr),
	modifier(mod)
{
}

double Exclusion::Bisect(double S, double E, double &n)
{
	double fS = Zero(S);
	double fE = Zero(E);

	if (fS * fE > 0)
	{
		n = fE + threshold;
		return E;
	}

	while (fabs(S - E) > 1e-6)
	{
		//std::cout << "Ranging " << S << " (" << fS << ") - "
		//			<< E << " (" << fE << ")" << std::endl;
		double M = (S + E) / 2.0;
		double fM = Zero(M);

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

	n = NumberEvents( (S + E) / 2.0 );
	return (S + E) / 2.0;
}

//the two intervals are [S : M] and [M : E]
bool Exclusion::FindInterval(double S, double &M, double E)
{
	double numS = Zero(S);
	double numE = Zero(E);

	for (int i = 0; i < 8; ++i)	//max 2^8 = 256 intervals
	{
		double split = (E - S) / pow (2, i+1);
		for (int n = 1; n < pow(2, i+1); n += 2)
		{
			double numM = Zero(S + n*split);

			if (numS * numM < 0)
			{
				M = S + n * split;
				return true;
			}
		}
	}

	return false;
}

double Exclusion::Zero(double lu2)
{
	return NumberEvents(lu2) - threshold;
}

double Exclusion::NumberEvents(double lu2)
{
	double ue = UeFlag ? pow(10.0, 0.5 * lu2) : 0.,
	       um = UmFlag ? pow(10.0, 0.5 * lu2) : 0.,
	       ut = UtFlag ? pow(10.0, 0.5 * lu2) : 0.;

	return modifier * engine->MakeSampler(box, horn, ue, um, ut);
}

void Exclusion::SetThreshold(double T)
{
	threshold = T;
}
