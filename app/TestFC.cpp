#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <getopt.h>

#include "tools.h"
#include "analysis.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"background", 	required_argument, 	0, 'b'},
		{"confidence", 	required_argument, 	0, 'C'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	int b = 0;
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "b:C:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'b':
				b = strtol(optarg, NULL, 10);
				break;
			case 'C':
				CL = strtod(optarg, NULL);
				CL /= 100.0;
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}

	int nA, nB;
	double s = Belt(b, CL, nA, nB);
	
	std::cout << "For " << CL*100.0 << " %, mean signal is " << s << " from " << nA << " to " << nB << std::endl;

	return 0;
}
