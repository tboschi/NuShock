#include "physics/Productions.h"

namespace Production {

	// PRODUCTION STUFF

	Channel fromString(std::string chan) {
		if (chan == "MuonE")
			return Channel::MuonE;
		else if (chan == "MuonM")
			return Channel::MuonM;
		else if (chan == "TauEE")
			return Channel::TauEE;
		else if (chan == "TauET")
			return Channel::TauET;
		else if (chan == "TauMM")
			return Channel::TauMM;
		else if (chan == "TauMT")
			return Channel::TauMT;
		else if (chan == "TauPI")
			return Channel::TauPI;
		else if (chan == "Tau2PI")
			return Channel::Tau2PI;
		else if (chan == "PionE")
			return Channel::PionE;
		else if (chan == "PionM")
			return Channel::PionM;
		else if (chan == "KaonE")
			return Channel::KaonE;
		else if (chan == "KaonM")
			return Channel::KaonM;
		else if (chan == "CharmE")
			return Channel::CharmE;
		else if (chan == "CharmM")
			return Channel::CharmM;
		else if (chan == "CharmT")
			return Channel::CharmT;
		else if (chan == "Kaon0E")
			return Channel::Kaon0E;
		else if (chan == "Kaon0M")
			return Channel::Kaon0M;
		else if (chan == "KaonCE")
			return Channel::KaonCE;
		else if (chan == "KaonCM")
			return Channel::KaonCM;
		else
			throw std::invalid_argument("Production channel " + chan
							+ " is unknwon");
	}

	std::string toString(Channel chan) {
		switch (chan) {
			case Channel::MuonE:
				return "MuonE";
			case Channel::MuonM:
				return "MuonM";
			case Channel::TauEE:
				return "TauEE";
			case Channel::TauET:
				return "TauET";
			case Channel::TauMM:
				return "TauMM";
			case Channel::TauMT:
				return "TauMT";
			case Channel::TauPI:
				return "TauPI";
			case Channel::Tau2PI:
				return "Tau2PI";
			case Channel::PionE:
				return "PionE";
			case Channel::PionM:
				return "PionM";
			case Channel::KaonE:
				return "KaonE";
			case Channel::KaonM:
				return "KaonM";
			case Channel::CharmE:
				return "CharmE";
			case Channel::CharmM:
				return "CharmM";
			case Channel::CharmT:
				return "CharmT";
			case Channel::Kaon0E:
				return "Kaon0E";
			case Channel::Kaon0M:
				return "Kaon0M";
			case Channel::KaonCE:
				return "KaonCE";
			case Channel::KaonCM:
				return "KaonCM";
			default:
				throw std::invalid_argument("Production channel "
						+ toString(chan) + " is unknown");
		}
	}

	std::vector<Channel> Channels() {
		return { Channel::MuonE,
			 Channel::MuonM,
			 Channel::TauEE,
			 Channel::TauET,
			 Channel::TauMM,
			 Channel::TauMT,
			 Channel::TauPI,
			 Channel::Tau2PI,
			 Channel::PionE,
			 Channel::PionM,
			 Channel::KaonE,
			 Channel::KaonM,
			 Channel::CharmE,
			 Channel::CharmM,
			 Channel::CharmT,
			 Channel::Kaon0E,
			 Channel::Kaon0M,
			 Channel::KaonCE,
			 Channel::KaonCM };
	}

	std::vector<double> Masses(Channel chan) {
		switch (chan) {
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
				throw std::invalid_argument("Production channel "
						+ toString(chan) + " is unknown");
		}
	}

	std::vector<int> Pdgs(Channel chan) {
		switch(chan) {
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
				throw std::invalid_argument("Production channel "
						+ toString(chan) + " is unknwon");
		}
	}

	unsigned int Ns(Channel chan) {
		switch(chan) {
			case Channel::MuonE:
				return 3;
			case Channel::MuonM:
				return 3;
			case Channel::TauEE:
				return 3;
			case Channel::TauET:
				return 3;
			case Channel::TauMM:
				return 3;
			case Channel::TauMT:
				return 3;
			case Channel::TauPI:
				return 2;
			case Channel::Tau2PI:
				return 3;
			case Channel::PionE:
				return 2;
			case Channel::PionM:
				return 2;
			case Channel::KaonE:
				return 2;
			case Channel::KaonM:
				return 2;
			case Channel::CharmE:
				return 2;
			case Channel::CharmM:
				return 2;
			case Channel::CharmT:
				return 2;
			case Channel::Kaon0E:
				return 3;
			case Channel::Kaon0M:
				return 3;
			case Channel::KaonCE:
				return 3;
			case Channel::KaonCM:
				return 3;
			default:
				throw std::invalid_argument("Production channel "
						+ toString(chan) + " is unknwon");
		}
	}
}
