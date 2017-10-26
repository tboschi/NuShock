#include "Nucleon.h"

Nucleon::Nucleon(bool Neutrino, bool Nucleon) :
	M_Proton(Const::fMProton),
	M_Neutron(Const::fMNeutron)
{
	Mh_prev = -1.0;
	El_prev = -1.0;

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
	if (!Event->SetDecay(ss, 2, Mass))
	{
		std::cout << Pl->E() << "\t";
		std::cout << "WARNING! Kinematics don't permit this process" << std::endl;
		std::cout << "Maximum heavy neutrino mass is " << ss.Mag() - Mn << " GeV" << std::endl;
	}

	GeneratePS();				//Set the 4vecs
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
			//( pow((s() - Ml*Ml - Mn*Mn), 2) - 4 * Ml*Ml * Mn*Mn );
	}
	else return 0.0;
}

double Nucleon::dSigmadOmega(int Neut)	//differential cross section (dS/dQ2) in CM frame	in GeV-4
{
	if (IsEnergyConserved())
	{
		return  Amp2(Neut) / 
			( 64.0 * Const::fPi2 * s() ) * 
			( pow((s() - Mh*Mh - Mn*Mn), 2) - 4 * Mh*Mh * Mn*Mn ) / 
			( pow((s() - Ml*Ml - Mn*Mn), 2) - 4 * Ml*Ml * Mn*Mn );
	}
	else return 0.0;
}

double Nucleon::SigmaTot(int Neut)	//differential cross section (dS/dQ2) in any frame	in GeV-2
{
	//std::cout << "Energy " << IsEnergyConserved() << std::endl;
	//std::cout << "s " << s() << "\tt " << t() << "\tu " << u() << "\tstu " << 2*Mn*Mn + Mh*Mh + Ml*Ml << std::endl;
	if (IsChanged())
	{
		double Q2min, Q2max;
		double dQ2 = Q2Lim(Q2min, Q2max);
		//std::cout << "integration " << Q2min << "\t" << Q2max << "\tdQ2 " << dQ2 << std::endl;
	
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
		if (A() < 0 || B() < 0 || C() < 0)
			std::cout << "help" << A() << "\t" << B() << "\t" << C() <<std::endl;
		return  8 * Const::fGF2 * pow(Mn, 4) *
			( A() +
			  x * B() * (s() - u()) / (Mn*Mn) +
			  C() * pow((s() - u()) / (Mn*Mn), 2) );
	}
	else return 0.0;
}

/*
//Maximum values of dS for MC purposes
double Nucleon::MaxSigma()
{
	double temx = x();
	double temy = y();

	if (IsChanged())
	{
		fMax = 0.0;

		//Phasespace coordinates are already checked by ddGamma
		//for (double ix = 0.0; ix <= 2.0; ix += 2.0/Kine::Loop)
		for (double ix = 0.0; ix <= 2.0; ix += 2.0/100)
		{
			SetX(ix);
			//for (double iy = 0.0; iy <= 2.0; iy += 2.0/Kine::Loop)
			for (double iy = 0.0; iy <= 2.0; iy += 2.0/100)
			{
				SetY(iy);
				double Gam = ddGamma();
				if (fMax < Gam)
				       fMax = Gam;	
			}
		}
	}

	SetX(temx);
	SetY(temy);

	return fMax;
}
*/

//for boole integeration with Kine
double Nucleon::Variable(double dt)
{
	double A, B;
	double dQ2 = Q2Lim(A, B);
	if (dQ2 > 0)
	{
		SetQ2(dQ2*dt + A);
		return dQ2*dSigmadQ2(1);
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

double Nucleon::Q2Lim(double &Min, double &Max)	//y integration limits
{
	double XX = (El+Mn)*(2*El*Mn + Ml*Ml + Mh*Mh);
	double AA = El*El - Ml*Ml;
	double BB = (Ml*Ml - Mh*Mh) * (Ml*Ml - Mh*Mh + 4*El*Mn) + 4*Mn*Mn*(El*El - Mh*Mh);
	double DD = 2*(2*El*Mn + Mn*Mn + Ml*Ml);

	//std::cout << "El " << El << "\tMl " << Ml << "\tMh " << Mh << "\t";
	//std::cout << "XX " << XX << "\tAA " << AA << "\tBB " << BB << "\tDD " << DD << std::endl;

	//safety of sqrt
	AA = AA < 0 ? 0.0 : AA;
	BB = BB < 0 ? 0.0 : BB;
	double Ehmin = (XX - sqrt(AA) * sqrt(BB)) / DD;
	double Ehmax = (XX + sqrt(AA) * sqrt(BB)) / DD;

	//std::cout << "Ehmin " << Ehmin << "\t Ehmax " << Ehmax << std::endl;
	Min = 2*Mn*(El - Ehmax);
	Max = 2*Mn*(El - Ehmin);

	//std::cout << "qmin " << Min << "\t qmax " << Max << std::endl;

	return Max - Min;
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
		fB = Q2() / (Mn*Mn) *
		     GA() * ( F1() + F2() );

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
	if (fF1 < 0)
		fF1 = (1.0/2.0 - Const::fSin2W) * 
		      (F1EM(1) - F1EM(0)) * Tau3() - 
		       Const::fSin2W * (F1EM(1) + F1EM(0)) -
		      F1S()/2.0;

	return pow(fF1, e);
}

double Nucleon::F2(int e)
{
	if (fF2 < 0)
		fF2 = (1.0/2.0 - Const::fSin2W) * 
		      (F2EM(1) - F2EM(0)) * Tau3() - 
		       Const::fSin2W * (F2EM(1) + F2EM(0)) -
		       F2S()/2.0;

	return pow(fF2, e);
}

double Nucleon::GA(int e)
{
	if (fGA < 0)
		fGA = (Tau3() * GAEM() - GAS())/ 2.0;

	return pow(fGA, e);
}

double Nucleon::F1EM(bool N, int e)
{
	if (fF1EM < 0)
	{
		double Mass = N ? M_Proton : M_Neutron;
		fF1EM = (GEsachs(N) + GMsachs(N) * Q2()/(4*Mass*Mass)) / 
			(1 + Q2()/(4*Mass*Mass));
	}
	
	return pow(fF1EM, e);
}

double Nucleon::F2EM(bool N, int e)
{
	if (fF2EM < 0)
	{
		double Mass = N ? M_Proton : M_Neutron;
		fF2EM = (GMsachs(N) + GEsachs(N)) / 
	        	(1 + Q2()/(4*Mass*Mass));
	}
	
	return pow(fF2EM, e);
}

double Nucleon::GAEM(int e)
{
	if (fGAEM < 0)
	{
		double eta = 0.12;	//GENIE parametrisation for correction

		double MA = Const::fMA;
		fGAEM = Const::fGA0 / pow(1+ Q2()/(MA*MA), 2) * (1 + eta);
	}

	return pow(fGAEM, e);
}

double Nucleon::F1S(int e)	//negligible atm
{
	if (fF1S < 0)
		fF1S = 0.0;

	return pow(fF1S, e);
}

double Nucleon::F2S(int e)	//negligible atm
{
	if (fF2S < 0)
		fF2S = 0.0;

	return pow(fF2S, e);
}

double Nucleon::GAS(int e)	//negligible atm
{
	if (fGAS < 0)
		fGAS = 0.0;

	return pow(fGAS, e);
}

double Nucleon::GEsachs(bool N, int e)
{
	if (fGEsachs < 0)
		fGEsachs = N ? GDipole() : 0;

	return pow(fGEsachs, e);
}

double Nucleon::GMsachs(bool N, int e)
{
	if (fGMsachs < 0)
		fGMsachs = N ? Const::fMagMuP * GDipole() : Const::fMagMuN * GDipole();

	return pow(fGMsachs, e);
}

double Nucleon::GDipole()
{
	if (fGDipole < 0)
	{
		double MV = Const::fMV;
		fGDipole = pow(1 + Q2()/(MV*MV), -2);
	}

	return fGDipole;
}

//Mandelstam variables
double Nucleon::s(double e)
{
	//TLorentzVector S(Pl->Vect()+Pi->Vect(), Pl->E()+Pi->E());
	//return S.M2();
	return pow(Ml*Ml + Mn*Mn + 2*El*Pi->E() - 2*Pl->Vect()*Pi->Vect(), e);
}

double Nucleon::t(double e)
{
	//TLorentzVector T(Pl->Vect()-Ph->Vect(), Pl->E()-Ph->E());
	//return T.M2();
	return pow(Ml*Ml + Mh*Mh - 2*El*Ph->E() + 2*Pl->Vect()*Ph->Vect(), e);
}

double Nucleon::u(double e)
{
	//TLorentzVector U(Pl->Vect()-Pf->Vect(), Pl->E()-Pf->E());
	//return U.M2();
	return pow(Ml*Ml + Mn*Mn - 2*El*Pf->E() + 2*Pl->Vect()*Pf->Vect(), e);
}

double Nucleon::Q2()
{
	return -t();
	//return q2;
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
void Nucleon::SetQ2(double X)
{
	double Eh = El - X/(2*Mn);

	if (Eh*Eh > Mh*Mh)	//Ph->Vect() is not zero, there is an angle
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

		TLorentzVector NewReco(Pl->Vect(), El + Mn - Eh);
		SetRecoil(NewReco);
	}

	ResetFormFactors();

	/*
	std::cout << "B: Q2 " << X << "\t";
	std::cout << "Eh " << Eh << "\t"; 
	std::cout << "Eh " << Ph->E() << "\t"; 
	std::cout << "m3 " << Ph->M() << "\t"; 
	std::cout << "cos " << cosT << "\t";
	std::cout << "thet " << Acos(cosT) << "\t";
	std::cout << "px " << Ph->Px() << "\tpy " << Ph->Py() << "\tpz " << Ph->Pz() << std::endl;
	*/
}

void Nucleon::SetNeutrino(bool B)
{
	Nu = -1 + 2*B;

	TLorentzVector nu(0, 0, 20, 20);
	SetProbe(nu);
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

bool Nucleon::IsEnergyConserved()
{
	if (fabs(s() + t() + u() - 2*Mn*Mn - Ml*Ml - Mh*Mh) < 1e-9)
		return true;
	else return false;
}

double Nucleon::Acos(double X)
{
	if (X > 1.0) X = 1.0;
	else if (X < -1.0) X = -1.0;
	return acos(X);
}

bool Nucleon::IsChanged()
{
	bool Ret = (fabs(Mh - Mh_prev) > 1e-9 ||
		    fabs(El - El_prev) > 1e-9);

	Mh_prev = Mh;
	El_prev = El;

	return Ret;
}

void Nucleon::ResetFormFactors()
{
	fA = -999.9;
	fB = -999.9;
	fC = -999.9;
	fF1 = -1.0;
	fF2 = -1.0;
	fGA = -1.0;
	fF1EM = -1.0;
	fF2EM = -1.0;
	fGAEM = -1.0;
	fF1S = -1.0;
	fF2S = -1.0;
	fGAS = -1.0;
	fGEsachs = -1.0;
	fGMsachs = -1.0;
	fGDipole = -1.0;
}
