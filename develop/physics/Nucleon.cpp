#include "Nucleon.h"

Nucleon::Nucleon(bool Neutrino, bool Nucleon) :
	M_Proton(Const::fMProton),
	M_Neutron(Const::fMNeutron)
{
	Mh_prev = -1.0;
	El_prev = -1.0;
	fQ2Min = -999.9;
	fQ2Max = -999.9;
	ResetFormFactors();

	Pl = new TLorentzVector();	//light neutrino
	Pi = new TLorentzVector();	//initial nucleon (rest)
	Ph = new TLorentzVector();	//heavy neutrino
	Pf = new TLorentzVector();	//final nucleon

	SetNeutrino(Neutrino);
	SetNucleon(Nucleon);

	Event = new TGenPhaseSpace;
}

void Nucleon::SetPS(double M_Sterile)		//PS setup for MC
{
	Mh = M_Sterile;
	double Mass[2] = {Mh, Mn};
	TLorentzVector ss = *Pl + *Pi;
	if (Event->SetDecay(ss, 2, Mass))
		GeneratePS();
}

double Nucleon::GeneratePS()	//generate event
{
	double Weight = Event->Generate();

	Ph->SetVect(Event->GetDecay(0)->Vect());
	Ph->SetE(Event->GetDecay(0)->E());

	Pf->SetVect(Event->GetDecay(1)->Vect());
	Pf->SetE(Event->GetDecay(1)->E());

	return Weight;
}

double Nucleon::dSigmadQ2(int Neut)	//differential cross section (dS/dQ2) in any frame
{
	if (IsEnergyConserved())
	{
		return  Amp2(Neut) / 
			( 64.0 * Const::fPi ) / 
			( pow(Pl->Dot(*Pi), 2) - Mn*Mn * Ml*Ml );
	}
	else return 0.0;
}

double Nucleon::dSigmadOmega(int Neut)	//differential cross section (dS/dQ2) in CM frame	in GeV-4
{
	if (IsEnergyConserved())
	{
		return  Amp2(Neut) / 
			( 64.0 * Const::fPi2 * s() ) * 
			( pow(Ph->Dot(*Pi), 2) - Mn*Mn * Mh*Mh ) /
			( pow(Pl->Dot(*Pi), 2) - Mn*Mn * Ml*Ml );
	}
	else return 0.0;
}

double Nucleon::SigmaTot(int Neut)	//differential cross section (dS/dQ2) in any frame	in GeV-2
{
	if (IsChanged())
	{
		double Q2min, Q2max;
		double dQ2 = Q2Lim(Q2min, Q2max);
	
		if (dQ2 > 0)
		{
			double (Nucleon::*dS)(int) = &Nucleon::dSigmadQ2;
			fSigmaTot = Integrate(dS, Q2min, Q2max, Neut);
		}
		else 
			fSigmaTot = 0.0;
	}

	return fSigmaTot;
}

double Nucleon::Amp2(int Neut)		//Unpolarised amplitude
{
	if (IsEnergyConserved())
	{
		int x = (Neut == 0) ? Sign() : Neut;

		return  8 * Const::fGF2 * pow(Mn, 4) *
			( A() +
			  x * B() * (s() - u()) / (Mn*Mn) +
			  C() * pow((s() - u()) / (Mn*Mn), 2) );
	}
	else return 0.0;
}

	
//form factors

double Nucleon::A()
{
	if (fA < -999.0)
	       fA = (Q2() + Ml*Ml + Mh*Mh) / (Mn*Mn) *
		    ((1 + Q2()/(4*Mn*Mn)) * GA(2) -
     		     (1 - Q2()/(4*Mn*Mn)) * (F1(2) - F2(2)* Q2()/(4*Mn*Mn)) +
     		     Q2()/(Mn*Mn) * F1()*F2()) -
     		    (Q2() * (Ml*Ml + Mh*Mh) + pow((Ml*Ml - Mh*Mh), 2)) / (4*pow(Mn,4)) *
     		    (F1(2) + F2(2) + 2*F1()*F2() + GA(2));

	return fA;
}

double Nucleon::B()
{
	if (fB < -999.0)
		fB = Q2() / (Mn*Mn) * GA() * ( F1() + F2() );

	return fB;
}

double Nucleon::C()
{
	if (fC < -999.0)
		fC = 1.0/4.0 * ( GA(2) + F1(2) + F2(2) * Q2()/(4*Mn*Mn) );

	return fC;
}

double Nucleon::F1(int e)
{
	if (fF1 < -999.0)
		fF1 = (1.0/2.0 - Const::fSin2W) * (F1EM(1) - F1EM(0)) * Tau3() - 
		      F1S()/2.0 - Const::fSin2W * (F1EM(1) + F1EM(0));

	return pow(fF1, e);
}

double Nucleon::F2(int e)
{
	if (fF2 < -999.0)
		fF2 = (1.0/2.0 - Const::fSin2W) * (F2EM(1) - F2EM(0)) * Tau3() - 
		      F2S()/2.0 - Const::fSin2W * (F2EM(1) + F2EM(0));

	return pow(fF2, e);
}

double Nucleon::GA(int e)
{
	if (fGA < -999.0)
		fGA = (Tau3() * GAEM() - GAS())/ 2.0;

	return pow(fGA, e);
}

double Nucleon::F1EM(bool N, int e)
{
	double Mass = N ? M_Proton : M_Neutron;
	fF1EM = (GEsachs(N) + GMsachs(N) * Q2()/(4*Mass*Mass)) / 
		(1 + Q2()/(4*Mass*Mass));
	
	return pow(fF1EM, e);
}

double Nucleon::F2EM(bool N, int e)
{
	double Mass = N ? M_Proton : M_Neutron;
	fF2EM = (GMsachs(N) + GEsachs(N)) / 
	       	(1 + Q2()/(4*Mass*Mass));
	
	return pow(fF2EM, e);
}

double Nucleon::GAEM(int e)
{
	double eta = 0.12;	//GENIE parametrisation for correction

	double MA = Const::fMA;
	fGAEM = Const::fGA0 / pow(1+ Q2()/(MA*MA), 2) * (1 + eta);

	return pow(fGAEM, e);
}

double Nucleon::F1S(int e)	//negligible atm
{
	fF1S = 0.0;

	return pow(fF1S, e);
}

double Nucleon::F2S(int e)	//negligible atm
{
	fF2S = 0.0;

	return pow(fF2S, e);
}

double Nucleon::GAS(int e)	//negligible atm
{
	fGAS = 0.0;

	return pow(fGAS, e);
}

double Nucleon::GEsachs(bool N, int e)
{
	fGEsachs = N ? GDipole() : 0;

	return pow(fGEsachs, e);
}

double Nucleon::GMsachs(bool N, int e)
{
	fGMsachs = N ? Const::fMagMuP * GDipole() : Const::fMagMuN * GDipole();

	return pow(fGMsachs, e);
}

double Nucleon::GDipole()
{
	if (fGDipole < -999.0)
	{
		double MV = Const::fMV;
		fGDipole = pow(1 + Q2()/(MV*MV), -2);
	}

	return fGDipole;
}

//Mandelstam variables
double Nucleon::s()	//1+2
{
	TLorentzVector S = *Pl + *Pi;
	return S.M2();
}

double Nucleon::s_()	//3+4
{
	TLorentzVector S = *Ph + *Pf;
	return S.M2();
}

double Nucleon::t()	//1-3
{
	TLorentzVector T = *Pl - *Ph;
	return T.M2();
}

double Nucleon::t_()	//4-2
{
	TLorentzVector T = *Pf - *Pi;
	return T.M2();
}

double Nucleon::u()	//1+4
{
	TLorentzVector U = *Pl - *Pf;
	return U.M2();
}

double Nucleon::u_()	//3-2
{
	TLorentzVector U = *Ph - *Pi;
	return U.M2();
}

double Nucleon::Q2()
{
	return -t();
}	

double Nucleon::Q2Lim(double &Min, double &Max)	//y integration limits
{
	if (fQ2Min < -999.0 || fQ2Max < -999.0 || IsChanged())
	{
		double XX = (El+Mn)*(2*El*Mn + Ml*Ml + Mh*Mh);
		double AA = El*El - Ml*Ml;
		double BB = (Ml*Ml - Mh*Mh) * (Ml*Ml - Mh*Mh + 4*El*Mn) + 4*Mn*Mn*(El*El - Mh*Mh);
		double DD = 2*(2*El*Mn + Mn*Mn + Ml*Ml);
	
		//safety of sqrt
		AA = AA < 0 ? 0.0 : AA;
		BB = BB < 0 ? 0.0 : BB;
		double Ehmin = (XX - sqrt(AA) * sqrt(BB)) / DD;
		double Ehmax = (XX + sqrt(AA) * sqrt(BB)) / DD;
	
		fQ2Min = 2*Mn*(El - Ehmax);
		fQ2Max = 2*Mn*(El - Ehmin);
	}

	if (IsAllowed())
	{
		Min = fQ2Min;
		Max = fQ2Max;
		return fQ2Max - fQ2Min;
	}
	else
		return -1.0;
}

void Nucleon::SetQ2(double X)
{
	double Eh = El - X/(2*Mn);

	if (Eh > Mh)	//Ph->Vect() is not zero, there is an angle
	{
		double cosT = (2*El*Eh - Ml*Ml - Mh*Mh - X) /
		      (2 * Pl->Vect().Mag() * sqrt(Eh*Eh - Mh*Mh) );
		SetSterileE(Eh);
		Ph->SetTheta(Acos(cosT));

		TLorentzVector NewReco(Pl->Vect()-Ph->Vect(), El + Mn - Eh);
		SetRecoil(NewReco);
	}
	else 	//Ph->Vect() is zero
	{
		SetSterileE(Mh);

		TLorentzVector NewReco(Pl->Vect(), El + Mn - Mh);
		SetRecoil(NewReco);
	}

	ResetFormFactors();
}

//Sign functions
int Nucleon::Tau3()	//+1 proton, -1 neutron
{
	return t3;
}

int Nucleon::Sign()	//T,+1 neutrino, F,-1 antineutrino
{
	return Nu;
}

double Nucleon::GetHeavyE()
{
	return Ph->E();
}	

//Setter
void Nucleon::SetNeutrino(bool B)
{
	TLorentzVector nu(0, 0, 20, 20);
	SetNeutrino(B, nu);
}

void Nucleon::SetNeutrino(bool B, TLorentzVector &nu)
{
	Nu = -1 + 2*B;

	SetProbe(nu);
}

void Nucleon::SetNucleon(bool B)
{
	t3 = -1 + 2*B;

	Mn = B ? M_Proton : M_Neutron;
	TLorentzVector A(0, 0, 0, Mn);
	SetTarget(A);
	SetRecoil(A);
}

void Nucleon::SetProbe(TLorentzVector &Probe)
{
	Pl->SetVect(Probe.Vect());
	Pl->SetE(Probe.E());

	Ml = Pl->M();
	El = Pl->E();
}

void Nucleon::SetProbeE(double dE)
{
	double dM = Pl->M();
	if (dE < dM)
	{
		Pl->SetE(dM);
		TVector3 Null;
		Pl->SetVect(Null);
	}
	else
	{
		Pl->SetE(dE);
		Pl->SetRho(sqrt(dE*dE - dM*dM));
	}
}

void Nucleon::SetProbeM(double dM)
{
	double dE = Pl->E();
	if (dE < dM)
	{
		Pl->SetE(dM);
		TVector3 Null;
		Pl->SetVect(Null);
	}
	else
		Pl->SetRho(sqrt(dE*dE - dM*dM));
}

void Nucleon::SetTarget(TLorentzVector &Targ)
{
	Pi->SetVect(Targ.Vect());
	Pi->SetE(Targ.E());
}

void Nucleon::SetTargetE(double dE)
{
	double dM = Mn;
	Pi->SetE(dE);
	if (dE < dM)
	{
		Pi->SetE(dM);
		TVector3 Null;
		Pi->SetVect(Null);
	}
	else
	{
		Pi->SetE(dE);
		Pi->SetRho(sqrt(dE*dE - dM*dM));
	}
}

void Nucleon::SetTargetM(double dM)
{
	double dE = Pi->E();
	if (dE < dM)
	{
		Pi->SetE(dM);
		TVector3 Null;
		Pi->SetVect(Null);
	}
	else
		Pi->SetRho(sqrt(dE*dE - dM*dM));
}

void Nucleon::SetSterile(TLorentzVector &Heavy)
{
	Ph->SetVect(Heavy.Vect());
	Ph->SetE(Heavy.E());
}

void Nucleon::SetSterileE(double dE)
{
	double dM = Ph->M();
	Ph->SetE(dE);
	if (dE < dM)
	{
		Ph->SetE(dM);
		TVector3 Null;
		Ph->SetVect(Null);
	}
	else
	{
		Ph->SetE(dE);
		Ph->SetRho(sqrt(dE*dE - dM*dM));
	}
}

void Nucleon::SetSterileM(double dM)
{
	
	double dE = Ph->E();
	if (dE < dM)
	{
		Ph->SetE(dM);
		TVector3 Null;
		Ph->SetVect(Null);
	}
	else
		Ph->SetRho(sqrt(dE*dE - dM*dM));
}

void Nucleon::SetRecoil(TLorentzVector &Reco)
{
	Pf->SetVect(Reco.Vect());
	Pf->SetE(Reco.E());
}

void Nucleon::SetRecoilE(double dE)
{
	double dM = Mn;
	Pf->SetE(dE);
	if (dE < dM)
	{
		Pf->SetE(dM);
		TVector3 Null;
		Pf->SetVect(Null);
	}
	else
	{
		Pf->SetE(dE);
		Pf->SetRho(sqrt(dE*dE - dM*dM));
	}
}

void Nucleon::SetRecoilM(double dM)
{
	double dE = Pf->E();
	if (dE < dM)
	{
		Pl->SetE(dE);
		TVector3 Null;
		Pf->SetVect(Null);
	}
	else
		Pf->SetRho(sqrt(dE*dE - dM*dM));
}

double Nucleon::Acos(double X)
{
	if (X > 1.0) X = 1.0;
	else if (X < -1.0) X = -1.0;
	return acos(X);
}

bool Nucleon::IsEnergyConserved()
{
	if (fabs(Pl->E() + Pi->E() - Ph->E() - Pf->E()) < 1e-12)
		return true;
	else return false;
}

bool Nucleon::IsAllowed()
{
	if (sqrt(s()) > Mh + Mn)
		return true;
	else return false;
}

bool Nucleon::IsChanged()
{
	bool Ret = (fabs(Mh - Mh_prev) > 1e-9 ||
		    fabs(El - El_prev) > 1e-9);

	Mh_prev = Mh;
	El_prev = El;

	if (Ret)
	{
		fQ2Min = -999.9;
		fQ2Max = -999.9;
	}

	return Ret;
}

void Nucleon::ResetFormFactors()
{
	fA = -999.9;
	fB = -999.9;
	fC = -999.9;

	fF1 = -999.9;
	fF2 = -999.9;
	fGA = -999.9;
	fF1EM = -999.9;
	fF2EM = -999.9;
	fGAEM = -999.9;
	fF1S = -999.9;
	fF2S = -999.9;
	fGAS = -999.9;
	fGEsachs = -999.9;
	fGMsachs = -999.9;
	fGDipole = -999.9;
}

//for boole integeration with Kine
double Nucleon::Variable(double dt)
{
	double A, B;
	double dQ2 = Q2Lim(A, B);
	if (dQ2 > 0)
	{
		SetQ2(dQ2*dt + A);
		return dQ2*dSigmadQ2(0);
	}
	else 
		return 0.0;
}

double Nucleon::Integrate(double (Nucleon::*FF)(int), double A, double B, int Neut)
{
	double a = 0, b = 0;
	double h = (B-A)/100;
	double Integral = 0;	//Boole's method for integration
	for (a = A; b + 1e-9 < B; a = b)
	{
		b = a + h;
		SetQ2(a);
		Integral += 7*(this->*FF)(Neut);

		SetQ2((3*a+b)/4.0);
		Integral += 32*(this->*FF)(Neut);

		SetQ2((a+b)/2.0);
		Integral += 12*(this->*FF)(Neut);

		SetQ2((a+3*b)/4.0);
		Integral += 32*(this->*FF)(Neut);

		SetQ2(b);
		Integral += 7*(this->*FF)(Neut);
	}	

	return Integral * h/90.0;
}

//for vegas integeration with Kine
void Nucleon::Integrand(const int *nDim, const double *x, const int *nComp, double f[])
{
	double A, B;
	double dQ2 = Q2Lim(A, B);
	if (dQ2 > 0)
	{
		if (A > B)
		{
			double tmp = B;
			B = A;
			A = tmp;
		}

		SetQ2(dQ2*x[0] + A);
		f[0] = dQ2*dSigmadQ2();
	}
	else 
		f[0] = 0.0;
}
