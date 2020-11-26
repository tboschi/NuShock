#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools.h"
#include "flux.h"
#include "detector.h"
#include "physics.h"
#include "analysis.h"

#include "omp.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"detconfig", 	required_argument,	0, 'd'},
		{"fluxconfig", 	required_argument,	0, 'f'},
		{"channel", 	required_argument,	0, 'c'},
		{"output", 	required_argument,	0, 'o'},
		{"threshold", 	required_argument,	0, 't'},
		{"massdepend", 	required_argument,	0, 'q'},
		{"efficiency", 	optional_argument,	0, 'w'},
		{"particle", 	no_argument,		0, 'P'},
		{"antipart", 	no_argument,		0, 'A'},
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	double CL = 0.90, scale;			//confidence level
	std::string output;

	while((iarg = getopt_long(argc,argv, "o:r:C:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'o':
				output.assign(optarg);
				break;
			case 'r':
				scale = std::strtod(optarg, NULL);
				break;
			case 'C':
				CL = std::strtod(optarg, NULL);
				break;
			case 'h':
				std::cout << "HELP" << std::endl;
				return 1;
			default:
				break;
		}
	}

	std::ofstream out(output.c_str());
	for (double LNC_dirac = 0; LNC_dirac < 1000; LNC_dirac+=5)
	{
		double LNV_dirac = scale * LNC_dirac;

		double LNC_major = LNC_dirac + LNV_dirac;
		double LNV_major = LNC_dirac + LNV_dirac;

		std::vector<int> vNp_dir, vNm_dir;
		std::vector<int> vNp_maj, vNm_maj;

		std::cout << "dirac " << LNC_dirac << std::endl;

		//out << "#";
		//out << LNC_dirac << "\t" << LNV_dirac << "\t"
		//    << LNC_major << "\t" << LNV_major << "\t";
		if (IsSensitiveToLNV(CL, vNp_dir, vNm_dir, vNp_maj, vNm_maj,
					LNC_dirac, LNV_dirac, LNC_major, LNV_major));
			//out << 1 << std::endl;
		//else
			//out << 0 << std::endl;

		out << LNC_dirac << "\t" << LNV_dirac;
		for (int i = 0; i < vNp_dir.size(); ++i)
			out << "\t" << vNp_dir[i];
		for (int i = 0; i < vNm_dir.size(); ++i)
			out << "\t" << vNm_dir[i];

		out << "\t" << LNC_major;
		for (int i = 0; i < vNp_maj.size(); ++i)
			out << "\t" << vNp_maj[i];
		for (int i = 0; i < vNm_maj.size(); ++i)
			out << "\t" << vNm_maj[i];
		out << std::endl;
	}

	/*
	for (double fhc = 0; fhc < 200; ++fhc)
		for (double rhc = 0; rhc < 200; ++rhc)
	{
		int imax;
		int d_fhc, d_rhc;
		double sum = ConfidenceRegion_sloppy(CL, imax, d_fhc, d_rhc,
							  fhc,  rhc);

		if (sum >= CL)
		out << fhc << "\t" << rhc << "\t" << d_fhc << "\t" << d_rhc << "\t" << sum << std::endl;
	}
	*/

	return 0;
}
