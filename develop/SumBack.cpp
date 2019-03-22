#include <iostream>
#include <cstdlib>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"

#include "Detector.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string ListFile, BaseName;
	unsigned int Bins = 10;
	std::string eFile, bFile, mFile, detConfig;
	
	while((iarg = getopt_long(argc,argv, "d:E:M:B:", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'E':
				eFile.assign(optarg);
				break;
			case 'M':
				mFile.assign(optarg);
				break;
			case 'B':
				bFile.assign(optarg);
				break;
			default:
				break;
		}
	
	}

	Detector * theBox = new Detector(detConfig);

	//double weight = theBox->Get("WeightLAr") + theBox->Get("WeightFGT");
	//double totPOT = theBox->Get("POT/s") * 1.0e7 * theBox->Get("Years") / 1.0e20;
	double weight = 1.0;
	double totPOT = 1.0;

	double totCC_E = weight * totPOT * 3.0e3;	//per total weight and total POT
	double totNC_E = weight * totPOT * 1.0e3;	//per total weight and total POT
	double tot_E = totCC_E + totNC_E;

	double totCC_M = weight * totPOT * 240.0e3;	//per total weight and total POT
	double totNC_M = weight * totPOT * 79.0e3;	//per total weight and total POT
	double tot_M = totCC_M + totNC_M;

	double totCC_B = weight * totPOT * 17.8e3;	//per total weight and total POT
	double totNC_B = weight * totPOT * 7.3e3;	//per total weight and total POT
	double tot_B = totCC_B + totNC_B;

	double ew = (3.0+1.0);
	double mw = (240.0+79.0);
	double bw = (17.8+7.3);
	double tot = ew+mw+bw;
	ew /= tot;
	mw /= tot;
	bw /= tot;

	double  back_m0300_m = 1.0e-6 * 1;
	double  back_m0300_e = 1.0e-6 * 0.3;
	double  back_m0300_b = 1.0e-6 * 9;

	double  back_m1000_m = 1.0e-6 * 24;
	double  back_m1000_e = 1.0e-6 * 0.3;
	double  back_m1000_b = 1.0e-6 * 146;

	double  back_e0150_m = 1.0e-6 * 0.3;
	double  back_e0150_e = 1.0e-6 * 0.3;
	double  back_e0150_b = 1.0e-6 * 0.3;

	double  back_e0300_m = 1.0e-6 * 0.3;
	double  back_e0300_e = 1.0e-6 * 4;
	double  back_e0300_b = 1.0e-6 * 0.3;

	double  back_e1000_m = 1.0e-6 * 0.3;
	double  back_e1000_e = 1.0e-6 * 101;
	double  back_e1000_b = 1.0e-6 * 0.3;

	std::cout << "Background for MPI in LNC\n";
	std::cout << "0.3Gev:\t" << back_m0300_m * totCC_M + back_m0300_e * totCC_E << "\t+";
	std::cout << "       \t" << back_m0300_m * totNC_M + back_m0300_e * totNC_E << "\t=";
	std::cout << "       \t" << back_m0300_m *   tot_M + back_m0300_e *   tot_E << std::endl;
	std::cout << "1.0Gev:\t" << back_m1000_m * totCC_M + back_m1000_e * totCC_E << "\t+";
	std::cout << "       \t" << back_m1000_m * totNC_M + back_m1000_e * totNC_E << "\t=";
	std::cout << "       \t" << back_m1000_m *   tot_M + back_m1000_e *   tot_E << std::endl;

	std::cout << "Background for MPI in LNV\n";
	std::cout << "0.3Gev:\t" << back_m0300_b * totCC_B  << "\t+";     //+ back_m0300_e * totCC_E << "\t+";
	std::cout << "       \t" << back_m0300_b * totNC_B  << "\t=";     //+ back_m0300_e * totNC_E << "\t=";
	std::cout << "       \t" << back_m0300_b *   tot_B  << std::endl; //+ back_m0300_e *   tot_E << std::endl;
	std::cout << "1.0Gev:\t" << back_m1000_b * totCC_B  << "\t+";     //+ back_m1000_e * totCC_E << "\t+";
	std::cout << "       \t" << back_m1000_b * totNC_B  << "\t=";     //+ back_m1000_e * totNC_E << "\t=";
	std::cout << "       \t" << back_m1000_b *   tot_B  << std::endl; //+ back_m1000_e *   tot_E << std::endl;

	std::cout << std::endl;

	std::cout << "Background for EPI in LNC\n";
	std::cout << ".15Gev:\t" << back_e0150_m * totCC_M + back_e0150_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e0150_m * totNC_M + back_e0150_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e0150_m *   tot_M + back_e0150_e *   tot_E << std::endl;
	std::cout << "0.3Gev:\t" << back_e0300_m * totCC_M + back_e0300_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e0300_m * totNC_M + back_e0300_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e0300_m *   tot_M + back_e0300_e *   tot_E << std::endl;
	std::cout << "1.0Gev:\t" << back_e1000_m * totCC_M + back_e1000_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e1000_m * totNC_M + back_e1000_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e1000_m *   tot_M + back_e1000_e *   tot_E << std::endl;

	std::cout << "Background for EPI in LNV\n";
	std::cout << ".15Gev:\t" << back_e0150_b * totCC_B  << "\t+";     //+ back_e0150_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e0150_b * totNC_B  << "\t=";     //+ back_e0150_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e0150_b *   tot_B  << std::endl; //+ back_e0150_e *   tot_E << std::endl;
	std::cout << "0.3Gev:\t" << back_e0300_b * totCC_B  << "\t+";     //+ back_e0300_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e0300_b * totNC_B  << "\t=";     //+ back_e0300_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e0300_b *   tot_B  << std::endl; //+ back_e0300_e *   tot_E << std::endl;
	std::cout << "1.0Gev:\t" << back_e1000_b * totCC_B  << "\t+";     //+ back_e1000_e * totCC_E << "\t+";
	std::cout << "       \t" << back_e1000_b * totNC_B  << "\t=";     //+ back_e1000_e * totNC_E << "\t=";
	std::cout << "       \t" << back_e1000_b *   tot_B  << std::endl; //+ back_e1000_e *   tot_E << std::endl;

	return 0;
}       
