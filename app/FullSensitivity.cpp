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

#include "tools/CardDealer.h"

class NumberEvents
{
	public:
		NumberEvents(Engine* tE, Detector *tD, std::vector<std::string> &chs,
				bool ue = false,
				bool um = false,
				bool ut = false,
				double thr = 0) :
			engine(tE),
			box(tD),
			channels(chs),
			Ue(ue),
			Um(um),
			Ut(ut),
			threshold(thr)
		{
		}

		double Number(double lu2)
		{
		}

		double operator()(double lu2, double t = -1) const
		{
			//std::cout << "computing number of events at " << lu2 << " :\t";
			double ue = Ue ? pow(10.0, 0.5 * lu2) : 0.;
			double um = Um ? pow(10.0, 0.5 * lu2) : 0.;
			double ut = Ut ? pow(10.0, 0.5 * lu2) : 0.;

			double tot = 0;
			for (int i = 0; i < channels.size(); ++i)
			{
				engine->SetDecay(channels[i]);
				tot += engine->MakeSampler(box, ue, um, ut);
			}
			//std::cout << tot << " vs " << threshold << std::endl;
			return tot - (t < 0 ? threshold : t);
		}

		void SetThreshold(double t)
		{
			if (t >= 0)
				threshold = t;
		}

	private:
		Engine *engine;
		Detector *box;
		std::vector<std::string> channels;
		double Ue, Um, Ut;
		double threshold;
};

void usage(char* argv0);
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
		{"efficiency", 	no_argument,		0, 'W'},
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
	
	//Initialize variables
	std::string detConfig, module, fluxConfig, outName;
	std::string channel = "ALL";
	std::string background;
	
	bool left = false, right = false;	//default unpolarised
	bool particle = false, antipart = false;	//default both for dirac, used in lnv studies
	bool dirac = false;				//default majorana neutrino
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	double CL = 0.90;			//confidence level

	while((iarg = getopt_long(argc,argv, "d:l:f:c:w:o:C:EMTLRrjAPh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'd':
				detConfig.assign(optarg);
				break;
			case 'l':
				module.assign(optarg);
				break;
			case 'f':
				fluxConfig.assign(optarg);
				break;
			case 'c':
				channel.assign(optarg);
				break;
			case 'w':
				background.assign(optarg);
				break;
			case 'o':
				outName.assign(optarg);		//outName is the output path
				break;
			case 'C':
				CL = std::strtod(optarg, NULL);
				break;
			case 'E':
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
				break;
			case 'L':
				left = true;
				right = false;
				break;
			case 'R':
				left = false;
				right = true;
				break;
			case 'r':
				dirac = true;
				break;
			case 'j':
				dirac = false;
				break;
			case 'P':
				particle = true;
				antipart = false;
				break;
			case 'A':
				particle = false;
				antipart = true;
				break;
			case 'h':
				usage(argv[0]);
				return 1;
			default:
				break;
		}
	}

	//constructing the neutrino, engine and detector classes

	//output file name
	//
	if (channel.empty())
	{
		std::cerr << "You have to select one channel using the -c flag" << std::endl;
		return 1;
	}

	if (!UeFlag && !UmFlag && !UtFlag)
	{
		std::cerr << "You have to select at least one mixing, -E -M or -T" << std::endl;
		return 1;
	}

	if (outName.find(".dat") == std::string::npos)
		outName += ".dat";

	outName.insert(outName.find(".dat"), channel);

	if (!background.empty())
		outName.insert(outName.find(".dat"), "_W");

	if (UeFlag)
		outName.insert(outName.find(".dat"), "_E");
	if (UmFlag)
		outName.insert(outName.find(".dat"), "_M");
	if (UtFlag)
		outName.insert(outName.find(".dat"), "_T");
	
	unsigned int opt0, optB, optHel;
	if (dirac)
	{
		opt0 = Neutrino::Dirac;
		optB = Neutrino::Dirac | Neutrino::Antiparticle;
		outName.insert(outName.find(".dat"), "_d");
	}
	else
	{
		opt0 = optB = Neutrino::Majorana;
		outName.insert(outName.find(".dat"), "_m");
	}

	if (left)
	{
		optHel = Neutrino::Left;	//-1 helicity
		outName.insert(outName.find(".dat"), "_L");
	}
	else if (right)
	{
		optHel = Neutrino::Right;	//+1 helicity
		outName.insert(outName.find(".dat"), "_R");
	}
	else
	{
		optHel = Neutrino::Unpolarised;
		outName.insert(outName.find(".dat"), "_U");
	}


	//neutrino classes and decay channels
	//
	Neutrino nu0(0, opt0 | optHel);
	Neutrino nuB(0, optB | optHel);
	nu0.SetDecayChannel(channel);
	nuB.SetDecayChannel(channel);
	double mass0 = nu0.DecayThreshold(channel);
	double massE = nu0.ProductionThreshold();
	if (mass0 > massE)	//limit masses, swap just because you never know I messed up somewhere else..
	{
		double tmp = massE;
		massE = mass0;
		mass0 = tmp;
	}

	std::vector<std::string> channels;
	if (channel == "ExpALL" && !background.empty())
	{
		if (UeFlag)
			channels.push_back("EPI");
		if (UmFlag)
			channels.push_back("MPI");
		channels.push_back("nPI0");
		channels.push_back("nEE");
		channels.push_back("nEM");
		channels.push_back("nMM");
	}
	else
		channels.push_back(channel);

	//Engine
	//
	//Engine::Current horn;
	//if (particle)			//for LNV
	//	horn = Engine::FHC;
	//else if (antipart)		//for LNV
	//	horn = Engine::RHC;
	//else
	//	horn = Engine::both;
	Engine *theEngine = new Engine(fluxConfig);
	theEngine->BindNeutrino("nu0", nu0, Engine::FHC);
	theEngine->BindNeutrino("nuB", nuB, Engine::RHC);

	//Detector class
	//
	if (!module.empty())
	{
		outName.insert(outName.find(".dat"), "_");
		outName.insert(outName.find(".dat"), module);
	}
	Detector *ndBox = new Detector(detConfig, module);

	//a collection of parabolas for each channel
	std::map<std::string, std::vector<double> > mCurve;
	std::map<std::string, std::vector<double> >::iterator ic;
	if (!background.empty())
	{
		std::string ferm = dirac ? "dirac" : "major";
		CardDealer threshold(background);
		for (int i = 0; i < channels.size(); ++i)
		{
			std::string key = channels[i] + "_";
			std::map<std::string, std::vector<double> > mBeta;
			if (!threshold.Get(key, mBeta))
				continue;

			if (module.empty())
			{
				ndBox->SetEfficiency(key + "LAr_" + ferm, channels[i]);
				ndBox->SetEfficiency(key + "FGT_" + ferm, channels[i]);
			}
			else
				ndBox->SetEfficiency(key + module + "_" + ferm, channels[i]);

			std::vector<double> alpha(3);
			for (ic = mBeta.begin(); ic != mBeta.end(); ++ic)
			{
				if (ic->first.find(ferm) == std::string::npos)
					continue;

				std::vector<double> beta = ic->second;
				for (int j = 0; j < alpha.size(); ++j)
					alpha[j] += beta[j];
			}
			mCurve[channels[i]] = alpha;
		}
	}
	//mcurve size is the same as channels size, if ever set

	std::ofstream out(outName.c_str());

	std::vector<double> masses, mixingBot, mixingTop, numberBot, numberTop;
	NumberEvents numEvts(theEngine, ndBox, channels, UeFlag, UmFlag, UtFlag);
	int grid = 500;
	double m0 = 0.01, mE = 2.0;
	double tb = -1, thr = 0;
	for (int i = 0; i < grid; ++i)
	{
		double mass = m0 * pow(mE/m0, i / double(grid));
		std::cout << "mass " << mass << std::endl;

		if (mass < mass0 || mass > massE)
			continue;

		theEngine->GetNeutrino("nu0").SetMass(mass);
		theEngine->GetNeutrino("nuB").SetMass(mass);

		nu0.SetMass(mass);
		nu0.SetMixings(UeFlag ? 1 : 0, UmFlag ? 1 : 0, UtFlag ? 1 : 0);
		double totalBranch = 0;
		if (mCurve.size())
			for (int c = 0; c < channels.size(); ++c)
				totalBranch += nu0.DecayBranch(channels[c]);

		double totalBack = 0;
		for (ic = mCurve.begin(); ic != mCurve.end(); ++ic) //map of channel - paramters
		{
			double nBack = nu0.DecayBranch(ic->first) / totalBranch *
				     (ic->second[0] + mass * (ic->second[1] + mass * ic->second[2]));
			totalBack += nBack < 0 ? 0 : nBack;
		}

		if (tb != totalBack)
		{
			thr = Belt(int(totalBack), CL);
			tb = totalBack;
		}

		theEngine->MakeFlux();
		theEngine->ScaleToDetector(ndBox);
		numEvts.SetThreshold(thr);

		////////////////////////////////////////

		double lU2bot = -16.0;
		double lU2top =   0.0;
		double lU2mid = findInterval(numEvts, lU2bot, lU2top);
		//std::cout << "range\t" << lU2bot << "\t" << lU2mid << "\t" << lU2top << std::endl;
		if (lU2bot < lU2mid && lU2mid < lU2top)
		{
			lU2bot = bisect(numEvts, lU2bot, lU2mid);
			lU2top = bisect(numEvts, lU2mid, lU2top);

			masses.push_back(mass);
			mixingBot.push_back(pow(10, lU2bot));
			mixingTop.push_back(pow(10, lU2top));
			numberBot.push_back(numEvts(lU2bot, 0));
			numberTop.push_back(numEvts(lU2top, 0));
		}
	}

	for (int i = 0; i < masses.size(); ++i)
		out << masses[i] << "\t" << mixingBot[i]
		    << "\t" << numberBot[i] << "\n";

	for (int i = masses.size(); i > 0; --i)
		out << masses[i-1] << "\t" << mixingTop[i-1]
		    << "\t" << numberTop[i-1] << "\n";

	if (masses.size())
		out << masses[0] << "\t" << mixingBot[0]
		    << "\t" << numberBot[0] << std::endl;

	out.close();

	return 0;
}
	
void usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -d,  --detconfig [CONFIG]" << std::endl;
	std::cout << "\t\tDetector configuration file" << std::endl;
	std::cout <<"\n  -f,  --fluxconfig [CONFIG]" << std::endl;
	std::cout << "\t\tFlux configuration file" << std::endl;
	std::cout <<"\n  -o,  --output [path/STRING]" << std::endl;
	std::cout << "\t\tBase name of output file" << std::endl;
	std::cout <<"\n  -c,  --channel [STRING]" << std::endl;
	std::cout << "\t\tDecay channel, defaul ALL" << std::endl;
	std::cout <<"\n  -t,  --threshold [FLOAT]" << std::endl;
	std::cout << "\t\tThreshold intercept. If -q not set, then is mass independent. Defaul 2.44" << std::endl;
	std::cout <<"\n  -q,  --massdepend [FLOAT]" << std::endl;
	std::cout << "\t\tMass dependance of threshold. Needs -t flag. Defaul 0" << std::endl;
	std::cout <<"\n  -E,  -M,  -T" << std::endl;
	std::cout << "\t\tSelect which mixing element. Multiple selection allowed" << std::endl;
	std::cout <<"\n  -L,  -R,  -U" << std::endl;
	std::cout << "\t\tSelect neutrino polarisation." << std::endl;
	std::cout <<"\n  --dirac,  --majorana" << std::endl;
	std::cout << "\t\tSelect fermionic nature of neutrino. Default Majorana" << std::endl;
	std::cout <<"\n  -P, -A {incompatible with --majorana}" << std::endl;
	std::cout << "\t\tSelect either neutrino or antineutrino components for dirac. Default, both are used." << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
