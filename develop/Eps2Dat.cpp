#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <getopt.h>

#include "Tools.h"

#include "TFile.h"
#include "TH1.h"

int main(int argc, char** argv)
{

	const struct option longopts[] = 
	{
		{"input", 	required_argument, 	0, 'i'},
		{"output", 	required_argument, 	0, 'o'},
		{"x1",	 	required_argument, 	0, '0'},
		{"x2",	 	required_argument, 	0, '1'},
		{"y1",	 	required_argument, 	0, '2'},
		{"y2",	 	required_argument, 	0, '3'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::ifstream InFile;
	std::ofstream OutFile;
	
	while((iarg = getopt_long(argc,argv, "i:o:0:1:2:3:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				InFile.open(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				std::cout << "Convert EPS file to ROOT" << std::endl;
				std::cout << "Usage : " << std::endl;
				std::cout << "Eps2Root [OPTIONS]" << std::endl;
				std::cout <<"\n  -i,  --input" << std::endl;
				std::cout << "\t\tInput file is a plain file" << std::endl;
				std::cout <<"\n  -r,  --root" << std::endl;
				std::cout << "\t\tOutput file is a ROOT file" << std::endl;
				std::cout <<"\n  -o,  --output" << std::endl;
				std::cout << "\t\tOutput file is a text file" << std::endl;
				std::cout <<"\n  -y,  --yaxis" << std::endl;
				std::cout << "\t\tDefine bottomline for y axix" << std::endl;
				std::cout <<"\n  -h,  --help" << std::endl;
				std::cout << "\t\tPrint this message and exit" << std::endl;
				return 1;
			default:
				return 1;
		}
	}
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::vector<double> TheX, TheY;

	std::string Line;
	std::stringstream ssL;
	double PosX = 0, PosY = 0, X, Y;
	double X1, X2, Y1, Y2;
	double x1, x2, y1, y2;

	while(std::getline(InFile, Line))
	{
		ssL.str("");
		ssL.clear();

		if (Line[0] == '#')
		{
			Line.erase(Line.begin());
			ssL << Line;
			ssL >> X1;
			ssL >> X2;
			ssL >> Y1;
			ssL >> Y2;
		}
		else if (Line[0] == '$')
		{
			Line.erase(Line.begin());
			ssL << Line;
			ssL >> x1;
			ssL >> x2;
			ssL >> y1;
			ssL >> y2;
		}
		else
		{
			ssL << Line;
			ssL >> X >> Y;
			PosX += X;
			PosY += Y;
			TheX.push_back(PosX);
			TheY.push_back(PosY);
		}
	}

	//x = a+b*X
	
	double bX = log(x2/x1)/(X2 - X1);
	double aX = log(x1) - X1 * bX;
	
	double bY = log(y2/y1)/(Y2 - Y1);
	double aY = log(y1) - Y1 * bY;

	for (int i = 0; i < TheX.size(); ++i)
	{
		TheX.at(i) = aX + bX*TheX.at(i);
		TheY.at(i) = aY + bY*TheY.at(i);
	}

	for (int i = 0; i < TheX.size(); ++i)
		Out << exp(TheX.at(i)) << "\t" << exp(TheY.at(i)) << std::endl;

	return 0;
}
