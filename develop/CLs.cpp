#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <getopt.h>

class Sorter
{
	private:
		std::vector<long double> vInd;
	public:
		Sorter(std::vector<long double>& vvect) : vInd(vvect) {}
		bool operator()(int i, int j) const { return vInd.at(i) < vInd.at(j); }
};

void Usage(char* argv0);
long double FactPower(long double P, unsigned int e);
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
				Usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	std::cout << "Mean background is " << b << std::endl;
	std::cout << "Confidence Level is " << 100*CL << "\%" << std::endl;

	double s = sqrt(b);			//from plot

	long double Sum_SB = 1.0L, Sum_B = 1.0L;
	while (1 - Sum_SB / Sum_B < CL)
	{
		s += 0.001;

		long double bl = b;
		long double sbl = bl + s;
		long double Ratio_Obs = FactPower(sbl, b) / FactPower(bl, b);

		Sum_SB = 0.0L, Sum_B = 0.0L;
		for (unsigned int n = 0; n < b+1; ++n)	//detected events
		{
			long double P_SB = FactPower(sbl, n);
			long double P_B = FactPower(bl, n);

			long double Ratio = P_SB / P_B;

			Sum_SB += P_SB;
			Sum_B += P_B;
		}
		//std::cout << "signal " << s << "\tsum " << sum << "\tA " << nA << "\tB " << nB << std::endl;
		//std::cout << s << "\tSB " << Sum_SB << "\tB " << Sum_B << "\t1-CLs " << 1-Sum_SB/Sum_B << "\tobs" << Ratio_Obs << std::endl;
	}

	std::cout << "Mean signal is " << s << std::endl;
	return 0;
}

long double FactPower(long double P, unsigned int e)
{
	long double Ret = exp(-P);
	while (e > 0)
	{
		Ret *= P;
		Ret /= e--;
	}
	return Ret;
}	

void Usage(char* argv0)
{
	std::cout << "Implementation of Feldman Cousins method" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -b,  --background" << std::endl;
	std::cout << "\t\tExpected mean background" << std::endl;
	std::cout <<"\n  -C,  --confidence" << std::endl;
	std::cout << "\t\tPercentage of confidence interval. Default 90\%, typical values are 95\% and 99\%" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
