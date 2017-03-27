#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string SourceName)
{
	SourceFile = new TFile(SourceName.c_str(), "OPEN");

	TotalFlux = (TH1D*) SourceFile->Get("htotalflux");
	hMuonPion = (TH1D*) SourceFile->Get("hmuonpiona");
	hMuonKaon = (TH1D*) SourceFile->Get("hmuonkaon");
	hElectronPion = (TH1D*) SourceFile->Get("helectronpion");
	hElectronKaon = (TH1D*) SourceFile->Get("helectronkaon");
	hElectronKaon3 = (TH1D*) SourceFile->Get("helectronkaon3");
	hMuonKaonOther = (TH1D*) SourceFile->Get("hmuonkaonother");
}

void FluxDriver::MakeSterileFlux(double M_Sterile)
{
	double Data[4];
	LoadFlux(Energy, Data);
	FindFlux();
	ci;
	TMatrixD VanderM(4,4);
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			VanderM[i][j] = pow(Data[i], j);
	TMatrixD Result = VanderM.Invert()*Data;

	hS_MuonPion = Um*Um* ShrockFactor(genie::constant::MPION,,genieMMU) * nuFlux.MuPi;
	sFlux.MuKa = Um*Um* the_shrock_factor(MKAON,ms,MMU) * nuFlux.MuKa;
	
	sFlux.EPi  = Ue*Ue*the_shrock_factor(MPION,ms,ME)*   nuFlux.EPi;
	sFlux.EKa  = Ue*Ue*the_shrock_factor(MKAON,ms,ME) * nuFlux.EKa;

	sFlux.MuKaOther = the_shrock_factor(MKAON,ms,MMU)* Um*Um*nuFlux.MuKaOther;
	
	sFlux.EKa3 = Ue*Ue*nuFlux.EKa3;

}

Flux FluxDriver::SampleEnergy()
{
		random generation;
	
}

double FluxDriver::ShrockFactor(double M_Meson, double M_Lepton, double M_Sterile) //Correct scaling of flux
{
	double dM_a = pow(M_Lepton/M_Meson, 2.0);
	double dM_i = pow(M_Sterile/M_Meson, 2.0);	

	if (M_Meson >= M_Lepton + M_Sterile)
		return ShrockRho(dM_a, dM_i)/(dM_a * pow(1-dM_a, 2.0));
	else return 0;
}

double FluxDriver::ShrockRho(double X double Y)
{
	return ShrockFM(X, Y)*sqrt(ShrockLambda(1, X, Y));
}

double FluxDriver::ShrockFM(double X, double Y)
{
	return X+Y - (X-Y)*(X-Y);
}

double FluxDriver::ShrockLambda(double X, double Y, double Z)
{
	return X*X + Y*Y + Z*Z - 2*(X*Y + X*Z + Y*Z);
}
