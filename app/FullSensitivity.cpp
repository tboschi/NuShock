#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <getopt.h>

#include "tools/CardDealer.h"
#include "tools/RootFinding.h"

#include "detector/Detector.h"
#include "detector/Driver.h"

#include "montecarlo/Sampler.h"

#include "physics/Const.h"
#include "physics/Neutrino.h"
#include "physics/DecayRate.h"
#include "physics/ProductionRate.h"

void usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"channel", 	required_argument,	0, 'c'},
		{"threshold", 	required_argument,	0, 't'},
		{"massdepend", 	required_argument,	0, 'q'},
		{"efficiency", 	no_argument,		0, 'W'},
		{"dirac", 	no_argument,		0, 'r'},
		{"majorana", 	no_argument,		0, 'j'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::string channel = "ALL";
	std::string background;
	
	std::string ferm;
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	bool effSet = false;

	while((iarg = getopt_long(argc,argv, "c:w:EMTrjh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'c':
				channel.assign(optarg);
				break;
			case 'w':
				effSet = true;
				if (optarg)
					background.assign(optarg);
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
			case 'r':
				ferm = "d";
				break;
			case 'j':
				ferm = "m";
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
		throw std::logic_error("FullSens: a channel must be selected with -c <channel>");

	if (!UeFlag && !UmFlag && !UtFlag)
		throw std::invalid_argument("FullSens: at least one mixing must be selected with -E -M or -T");

	if (ferm.empty())
		throw std::logic_error("FullSens: fermion type must be selected with --dirac or --majorana");

	CardDealer cd(argv[optind]);

	std::string config;
	if (!cd.Get("detector_configuration", config))
		throw std::logic_error("There is no detector configuration\n");
	Detector box(config);

	if (!cd.Get("flux_configuration", config))
		throw std::logic_error("There is no flux configuration\n");
	Driver drive(config);

	std::string outname;
	if (!cd.Get("output", outname))
		throw std::logic_error("FullSens: no output \".dat\" file specified in card");
	if (outname.find(".dat") == std::string::npos)
		outname += ".dat";

	outname.insert(outname.find(".dat"), channel);

	if (!background.empty())
		outname.insert(outname.find(".dat"), "_W");

	if (UeFlag)
		outname.insert(outname.find(".dat"), "_E");
	if (UmFlag)
		outname.insert(outname.find(".dat"), "_M");
	if (UtFlag)
		outname.insert(outname.find(".dat"), "_T");
	
	outname.insert(outname.find(".dat"), "_" + ferm);

	/*
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
	*/


	/*
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
	*/
	//mcurve size is the same as channels size, if ever set

	// from card
	size_t grid;	// from card
	if (!cd.Get("discretization", grid))
		grid = 500;
	std::vector<double> mass_range;
	if (!cd.Get("mass_range", mass_range))
		mass_range = {0.01, 2.00};

	std::vector<double> mix_range;
	if (!cd.Get("lmix_range", mix_range))
		mass_range = {-16., 0.};

	double tb = -1, thr = 2.44;

	Mixing mix(UeFlag, UmFlag, UtFlag);
	Decay::Channel chan = Decay::fromString(channel);

	// specific for channel and fermion type
	//if (!effSet)
		//box.LoadEfficiency(chan, ferm);

	// store first result
	bool first = true;
	double m0, u0;

	std::vector<std::pair<double, double> > up_limit;
	up_limit.reserve(grid);

	//Sampler sample(box, drive); //, Neutrino(mass, ferm), mix);
	std::ofstream out(outname.c_str());
	for (size_t g = 0; g < grid; ++g) {
		double mass = std::pow(mass_range[1] / mass_range[0],
				       g / double(grid)) * mass_range[0];
		//if (mass < mass_thr)
		//	continue;
		
		// decay to channel not allowed for this mass
		if (!DecayRate::IsAllowed(mass, chan))
			continue;

		std::cout << "mass " << mass << "\n";

		/*
		std::vector<Spectrum> spectra;
		spectra.reserve(fermtypes.size());
		for (const auto & t : fermtypes) {
			std::cout << "making\n";
			auto spec = drive.MakeSpectrum(Neutrino(mass, t), mix);
			if (!spec(mix))	// empty
				continue;
			spectra.push_back(std::move(spec));
		}
	
		// no spectrum available for this mass
		if (!spectra.size())
			continue;
		*/

		// store here all neutrinos needed for the analysis
		std::vector<Neutrino> nus;
		if (ferm == "m")
			nus = {Neutrino(mass, Neutrino::majorana)};
		else if (ferm == "d")
			nus = {Neutrino(mass, Neutrino::dirac | Neutrino::particle),
			       Neutrino(mass, Neutrino::dirac | Neutrino::antiparticle)};

		// create a function that generate a sample for given mixing
		auto sample = Sampler::Build(box, drive, chan, mix, std::move(nus));

		// create function that computes number of events using sampler
		// lx is mixing in log scale
		auto func = [&](double lx) -> double {
			auto hist = sample(mix * std::pow(10., lx / 2.));
			return (hist ? hist->Integral() : 0.) - thr;
		};

		/*
		for (size_t u = 0; u < grid; ++u) {
			double lmix = mix_range[0] + u * (mix_range[1] - mix_range[0]) / grid;
			std::cout << "mass " << mass << " vs " << std::pow(10., lmix) << "\n";
			out << mass << "\t" << std::pow(10., lmix) << "\t" << func(lmix) << "\n";
} */

		double half = RootFinding::CheckInterval(func, mix_range[0], mix_range[1], 1.e-4);
		double lim = RootFinding::BinarySearch(func, mix_range[0], half, 1.e-4);
		if (half > mix_range[0]) {
			out << mass << "\t" << std::pow(10., lim) << "\n";
			if (first) {
				m0 = mass;
				u0 = std::pow(10., lim);
				first = false;
			}
			if (half < mix_range[1]) {
				lim = RootFinding::BinarySearch(func, half, mix_range[1], 1.e-4);
				up_limit.emplace_back(mass, pow(10., lim));
			}
		}
	}

	for (auto i = up_limit.rbegin(); i != up_limit.rend(); ++i)
		out << i->first << "\t" << i->second << "\n";
	out << m0 << "\t" << up_limit.front().second << "\n";
	out << m0 << "\t" << u0 << "\n";

	/*
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
		*/

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
