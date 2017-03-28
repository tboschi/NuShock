#include "EventGenerator.h"

EventGenerator::EventGenerator(double MSterile, double Ue, double Um, double Ut)
{
	TheGamma = new DecayRates(MSterile, Ue, Um, Ut);
}

double EventGenerator::SetEnergy(double X)
{
	E_Sterile = X;
}

double EventGenerator::Probability()
{
	double Ratio = TheGamma->Branch("channel"); 
	double Total = TheGamma->Total();
	double Lorentz = M_Sterile/sqrt(E_Sterile*E_Sterile - M_Sterile*M_Sterile);
	return exp(-Total * Length / Lorentz) * (1-exp(- Total * Lambda / Lorentz) * Ratio;
}
