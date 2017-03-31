#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <sys/stat.h>

#include "EventGenerator.h"
#include "TRandom3.h"

void Usage();
inline bool Good (const char *File);

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"smconfig", 	required_argument,	0, 's'},
		{"flux", 	required_argument, 	0, 'f'},
		{"detector", 	required_argument, 	0, 'd'},
		{"output", 	required_argument, 	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ofstream OutFile;
	
//Initialize variables
	std::string SMConfig, FluxConfig, DetectorConfig;

	while((iarg = getopt_long(argc,argv, "s:f:d:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 's':
				std::cout << "us" << std::endl;
				if (Good(optarg))
				{
					SMConfig.assign(optarg);
					break;
				}
				else return 1;
			case 'f':
				std::cout << "uf" << std::endl;
				if (Good(optarg))
				{
					FluxConfig.assign(optarg);
					break;
				}
				break;
			case 'd':
				std::cout << "ud" << std::endl;
				if (Good(optarg))
				{
					DetectorConfig.assign(optarg);
					break;
				}
				break;
			case 'o':
				std::cout << "uo" << std::endl;
				OutFile.open(optarg);
				break;
			case 'h':
				std::cout << "uh" << std::endl;
				Usage();
				return 1;
			default:
				std::cerr << "Uknown parameter" << std::endl; 
				std::cerr << "Type Exclusion --help or -h to see usage" << std::endl; 
				return 1;
		}
	
	}

	//To have multiple output, handled by usage
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	//Create object from classes
	EventGenerator *gEvGen = new EventGenerator(SMConfig, FluxConfig, DetectorConfig);
	TRandom3 *RanGen = new TRandom3(0);

	int dim = 100;
	double ** Matrix = new double*[dim];
	for (int i = 0; i < dim; ++i)
		Matrix[i] = new double[dim];

	//Main body
	Out << "#Mass\tCoupling\tValue" << std::endl;
	double Ui;
	double Mj;
	for (int i = 0; i < dim; ++i)
	{
		Ui = i*gEvGen->GetUe()/dim;
		for (int j = 0; j < dim; ++j)
		{ 
			Mj = j*gEvGen->GetMSterile()/dim;
			Out << Mj << "\t" << Ui << "\t" << RanGen->Gaus() << std::endl;
		}
	}

	//Garbage collection
	delete gEvGen;

	return 0;
}
	
void Usage()
{
	std::cout << "Generate exclusion plot" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << "Exclusion [OPTIONS]" << std::endl;
	std::cout <<"\n  -s,  --smconfig" << std::endl;
	std::cout << "\t\tConfiguration file for the physical parameters. NO DEFAULT" << std::endl;
	std::cout <<"\n  -f,  --flux" << std::endl;
	std::cout << "\t\tFlux input ROOT file. NO DEFAULT" << std::endl;
	std::cout <<"\n  -d,  --detector" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\t(ROOT) output file for saving exclusion 2D plot" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and quit" << std::endl;
}

inline bool Good (const char *File)
{
	std::ifstream Infile(File);
	return Infile.good();
}
