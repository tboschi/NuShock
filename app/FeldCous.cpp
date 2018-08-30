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
		bool operator()(int i, int j) const { return vInd.at(i) > vInd.at(j); }
};

void Usage(char* argv0);
long double Poisson(long double s, unsigned int b, unsigned int e);
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

	std::vector<long double> vRatio, vP0;
	std::vector<unsigned int> vN, vI;

	int nA = 0, nB = 0;

	std::cout << "Mean background is " << b << std::endl;
	std::cout << "Confidence Level is " << 100*CL << "\%" << std::endl;

	/* WARNING: BLACK MAGIC USED HERE!!! */
	double s = sqrt(2.5*b);			//from plot
	//double lim = 70*log(b + 270) - 350;	//from plot
	double lim = 70*log(b + 270);	//from plot
	/* END OF BLACK MAGIC */

	int Width = int(b + lim - std::max(0, b-1));
	std::cout << "Belt limited from " << std::max(0, b-1) << " to " << int(b + lim) << " (width " << Width << ")" << std::endl;
	while (nA < b+1)
	{
		//std::cout << nB-nA << std::endl;
		s += 0.005;
		vRatio.clear();
		vP0.clear();
		vN.clear();
		vI.clear();

		unsigned int j = 0;
		//for (unsigned int n = std::max(0, int(s+b-10)); n < s+b+10; ++n, ++j)
		//std::cout << std::endl;
		//std::cout << "min " << std::max(0,b-3) << "\tmax " << s+b+50 << "\t" << vI.size() << std::endl;
		for (unsigned int n = std::max(0,b-1); n < b + lim; ++n, ++j)
		//for (unsigned int n = 140; n < 160; ++n, ++j)
		{
			long double P0 = Poisson(s, b, n);
			long double PA = Poisson(std::max(0, int(n-b)), b, n);

			long double Ratio = P0 / PA;

			vRatio.push_back(Ratio);
			vP0.push_back(P0);
			vN.push_back(n);
			vI.push_back(j);

			//std::cout << "n " << n << "\ts " << s << "\tu " << std::max(0,int(n-b));
			//std::cout << "\tP0 " << P0 << "\tPA " << PA;
			//std::cout << "\tR " << Ratio << std::endl;
		}
		std::sort(vI.begin(), vI.end(), Sorter(vRatio));

		double sum = 0;
		//std::cout << "size " << vN.size() << std::endl;
		nA = int(b + Width);
		nB = 0;
		for (int i = 0; i < vI.size(); ++i)
		{
			sum += vP0.at(vI.at(i));
		//std::cout << vI.at(i) << "\tn " << vN.at(vI.at(i)) << "\tR " << vRatio.at(vI.at(i)) << "\t" << vP0.at(vI.at(i)) << "\tsum " << sum << std::endl;

			//std::cout << vN.at(vI.at(i)) << "\t";
			if (nA > vN.at(vI.at(i)))
				nA = vN.at(vI.at(i));
			if (nB < vN.at(vI.at(i)))
				nB = vN.at(vI.at(i));

			if (sum > CL)
				break;
		}
		//std::cout << "signal " << s << "\tsum " << sum << "\tA " << nA << "\tB " << nB << "\tfrom " << std::max(0,b-5) << "\tto " << b+lim+20 << std::endl;
		//std::cout << s << "\tA " << nA << "\tB " << nB << std::endl;
	}
	std::cout << "Mean signal is " << s << ", with belt from " << nA-1 << " to " << nB << " (width " << nB - nA + 1 << ")" << std::endl;
	if (nB-nA+1 > Width)
		std::cout << "There is a belt issue!" << std::endl;
	return 0;
}

long double Poisson(long double s, unsigned int b, unsigned int e)
{
	unsigned int c = 0;
	long double P = exp(-s);
	while (e > 0)
	{
		if (c < b)
		{
			P *= exp(-1.0L);
			c++;
		}
		P *= (b+s);
		P /= e--;
	}
	return P;
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
