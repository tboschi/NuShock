#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <getopt.h>

void Usage(char* Name) {};
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"list", 	required_argument, 	0, 'l'},
		{"output", 	required_argument, 	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	std::string ListFile, BaseName;
	unsigned int Bins = 10;
	
	while((iarg = getopt_long(argc,argv, "l:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'l':
				ListFile.assign(optarg);
				break;
			case 'o':
				BaseName.assign(optarg);
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	
	}

	std::ifstream List(ListFile.c_str());

	double Size;
	std::string RootFile;
	std::vector<double> vSize;
	std::vector<std::string> vFile;
	std::vector<double> vVal{319, 4, 25.1};
	double W_FGT = 8.0, W_LAr = 50.0; 
	std::vector<unsigned int> vFGT_A, vFGT_B, vLAr_A, vLAr_B;
	std::vector<unsigned int> *vRef_A, *vRef_B;

	unsigned int A, B;
	std::string Line;
	std::stringstream ssL;
	while (std::getline(List, Line))
	{
		if (Line[0] == '#') continue;

		ssL.clear();
		ssL.str("");

		if (Line.find("LAr") != std::string::npos)
		{
			vRef_A = &vLAr_A;
			vRef_B = &vLAr_B;
			continue;
		}
		if (Line.find("FGT") != std::string::npos)
		{
			vRef_A = &vFGT_A;
			vRef_B = &vFGT_B;
			continue;
		}

		ssL << Line;
		ssL >> A >> B;

		vRef_A->push_back(A);
		vRef_B->push_back(B);
	}

	List.close();

	std::cout << "per component" << std::endl;
	for (unsigned int i = 0; i < vVal.size(); ++i)
	{
		std::cout << i << "\t" << double(W_FGT*vFGT_A.at(i) + W_LAr*vLAr_A.at(i))/(W_FGT+W_LAr) << "\t->";
		std::cout <<      "\t" << double(W_FGT*vFGT_B.at(i) + W_LAr*vLAr_B.at(i))/(W_FGT+W_LAr) << std::endl;
	}

	std::cout << "average" << "\t";
	double Avg_A = 0.0, Avg_B = 0.0, weight = 0.0;
	for (unsigned int i = 0; i < vVal.size(); ++i)
	{
		Avg_A += (W_FGT*vFGT_A.at(i)+W_LAr*vLAr_A.at(i))*vVal.at(i);
		Avg_B += (W_FGT*vFGT_B.at(i)+W_LAr*vLAr_B.at(i))*vVal.at(i);
		weight += vVal.at(i);
	}

	std::cout << Avg_A/weight/(W_FGT+W_LAr) << "\t->\t" << Avg_B/weight/(W_FGT+W_LAr) << std::endl;

	return 0;
}

