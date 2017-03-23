#include <FluxDriver.h>

FluxDriver::FluxDriver(std::string SourceName)
{
	fopen;
	vector;
}

FluxDriver::MakeFlux(double Energy)
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
}

FluxDriver::SampleEnergy()
{
	random generation;
	
}
