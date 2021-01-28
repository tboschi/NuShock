#include "montecarlo/Process.h"

/*
 * functions return vectors of masses or pdg codes for
 * processes of interest.
 * The order of particles for decay channels is in increasing mass
 * as it is the code name for the channel
 * The order of particles for production channels is also
 * in increasing mass, but the first particle is the parent particle
 */

namespace Process {
	std::vector<double> Mass(Channel::Name chan) {
		switch(chan) {
			//DECAYS
			case Channel::ALL:
			case Channel::nnn:
				return {Const::MNeutrino, Const::MNeutrino, Const::MNeutrino};
			case Channel::nGAMMA:
				return {Const::MNeutrino, Const::MPhoton};
			case Channel::ExpALL:
			case Channel::nEE:
				return {Const::MNeutrino, Const::MElectron, Const::MElectron};
			case Channel::nMM:
				return {Const::MNeutrino, Const::MMuon, Const::MMuon};
			case Channel::nEM:
				return {Const::MNeutrino, Const::MElectron, Const::MMuon};
			case Channel::nET:
				return {Const::MNeutrino, Const::MElectron, Const::MTau};
			case Channel::nMT:
				return {Const::MNeutrino, Const::MMuon, Const::MTau};
			case Channel::nPI0:
				return {Const::MNeutrino, Const::MPion0};
			case Channel::nETA:
				return {Const::MNeutrino, Const::MEta};
			case Channel::nETAi:
				return {Const::MNeutrino, Const::MEtai};
			case Channel::EPI:
				return {Const::MElectron, Const::MPion};
			case Channel::MPI:
				return {Const::MMuon, Const::MPion};
			case Channel::TPI:
				return {Const::MTau, Const::MPion};
			case Channel::EKA:
				return {Const::MElectron, Const::MKaon};
			case Channel::MKA:
				return {Const::MMuon, Const::MKaon};
			case Channel::ECHARM:
				return {Const::MElectron, Const::MD};
			case Channel::nRHO0:
				return {Const::MNeutrino, Const::MRho0};
			case Channel::nOMEGA:
				return {Const::MNeutrino, Const::MOmega};
			case Channel::nPHI:
				return {Const::MNeutrino, Const::MPhi};
			case Channel::ERHO:
				return {Const::MElectron, Const::MRho};
			case Channel::MRHO:
				return {Const::MMuon, Const::MRho0};
			case Channel::EKAx:
				return {Const::MElectron, Const::MKaonx};
			case Channel::MKAx:
				return {Const::MMuon, Const::MKaonx
				//PRODUCTION
			case Channel::MuonE:
				return {Const::MMuon, Const::MElectron, Const::MNeutrino};
			case Channel::MuonM:
				return {Const::MMuon, Const::MElectron, Const::MNeutrino};
			case Channel::TauEE:
				return {Const::MTau, Const::MElectron, Const::MNeutrino};
			case Channel::TauET:
				return {Const::MTau, Const::MElectron, Const::MNeutrino};
			case Channel::TauMM:
				return {Const::MTau, Const::MMuon, Const::MNeutrino};
			case Channel::TauMT:
				return {Const::MTau, Const::MMuon, Const::MNeutrino};
			case Channel::TauPI:
				return {Const::MTau, Const::MPion};
			case Channel::Tau2PI:
				return {Const::MTau, Const::MPion, Const::MPion0};
			case Channel::PionE:
				return {Const::MPion, Const::MElectron};
			case Channel::PionM:
				return {Const::MPion, Const::MMuon};
			case Channel::KaonE:
				return {Const::MKaon, Const::MElectron};
			case Channel::KaonM:
				return {Const::MKaon, Const::MMuon};
			case Channel::CharmE:
				return {Const::MDs, Const::MElectron};
			case Channel::CharmM:
				return {Const::MDs, Const::MMuon};
			case Channel::CharmT:
				return {Const::MDs, Const::MTau};
			case Channel::Kaon0E:
				return {Const::MKaon0, Const::MPion, Const::MElectron};
			case Channel::Kaon0M:
				return {Const::MKaon0, Const::MPion, Const::MMuon};
			case Channel::KaonCE:
				return {Const::MKaon, Const::MPion0, Const::MElectron};
			case Channel::KaonCM:
				return {Const::MKaon, Const::MPion0, Const::MMuon};
			default:
				throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
		}
	}

	std::vector<int> Pdg(Channel::Name chan) {
		switch(chan) {
			//DECAYS
			case Channel::ALL:
			case Channel::nnn:
				return {12, -12, 12};
			case Channel::nGAMMA:
				return {12, 22};
			case Channel::ExpALL:
			case Channel::nEE:
				return {12, -11, 11};
			case Channel::nMM:
				return {12, -13, 13};
			case Channel::nEM:
				return {12, -11, 13};
			case Channel::nET:
				return {12, -11, 15};
			case Channel::nMT:
				return {12, -13, 15};
			case Channel::nPI0:
				return {12, 111};
			case Channel::nETA:
				return {12, 221};
			case Channel::nETAi:
				return {12, 331};
			case Channel::EPI:
				return {11, 211};
			case Channel::MPI:
				return {13, 211};
			case Channel::TPI:
				return {15, 211};
			case Channel::EKA:
				return {11, 321};
			case Channel::MKA:
				return {13, 321};
			case Channel::ECHARM:
				return {12, 411};
			case Channel::nRHO0:
				return {12, 113};
			case Channel::nOMEGA:
				return {12, 223};
			case Channel::nPHI:
				return {12, 333};
			case Channel::ERHO:
				return {11, 213};
			case Channel::MRHO:
				return {13, 213};
			case Channel::EKAx:
				return {11, 9000321};
			case Channel::MKAx:
				return {13, 9000321};
			// PRODUCTION
			case Channel::MuonE:
				return {13, 11, 14};
			case Channel::MuonM:
				return {13, 11, -12};
			case Channel::TauEE:
				return {15, 11, 16};
			case Channel::TauET:
				return {15, 11, -12};
			case Channel::TauMM:
				return {15, 13, 16};
			case Channel::TauMT:
				return {15, 13, -14};
			case Channel::TauPI:
				return {15, 211};
			case Channel::Tau2PI:
				return {15, -211, 111};
			case Channel::PionE:
				return {211, 11};
			case Channel::PionM:
				return {211, 13};
			case Channel::KaonE:
				return {321, 11};
			case Channel::KaonM:
				return {321, 13};
			case Channel::CharmE:
				return {431, 11};
			case Channel::CharmM:
				return {431, 13};
			case Channel::CharmT:
				return {431, 15};
			case Channel::Kaon0E:
				return {130, 211, -11};
			case Channel::Kaon0M:
				return {130, 211, -13};
			case Channel::KaonCE:
				return {321, 111, -11};
			case Channel::KaonCM:
				return {321, 111, -13};
			default:
				throw std::invalid_argument("Channel " + Channel::toString(chan) + " unknown");
		}
	}

	std::vector<std::pair<double, int> > MassPdg(Channel::Name chan) {
		std::vector<double> mass = Mass(chan);
		std::vector<int> pdgs = Pdg(chan);
		if (mass.size() != pdgs.size())
			throw std::logic_error("Process: masses and pdgs do not match\n");

		std::vector<std::pair<double, int> > masspdg;
		masspdg.reserve(mass.size());
		std::transform(mass.begin(), mass.end(), pdgs.begin(),
				std::back_inserter(masspdg),
				[](double m, int p) { return std::make_pair(m, p); } );

		return masspdg;
	}
}
