#include "physics/Decays.h"

namespace Decay {

	// DECAY STUFF

	Channel fromString(std::string chan) {
		//if (chan == "ALL")
			//return Channel::ALL;
		if (chan == "nnn")
			return Channel::nnn;
		else if (chan == "nGAMMA")
			return Channel::nGAMMA;
		else if (chan == "nEE")
			return Channel::nEE;
		else if (chan == "nEM")
			return Channel::nEM;
		else if (chan == "nMM")
			return Channel::nMM;
		else if (chan == "nET")
			return Channel::nET;
		else if (chan == "nMT")
			return Channel::nMT;
		else if (chan == "nPI0")
			return Channel::nPI0;
		else if (chan == "EPI")
			return Channel::EPI;
		else if (chan == "MPI")
			return Channel::MPI;
		else if (chan == "TPI")
			return Channel::TPI;
		else if (chan == "EKA")
			return Channel::EKA;
		else if (chan == "MKA")
			return Channel::MKA;
		else if (chan == "nRHO0")
			return Channel::nRHO0;
		else if (chan == "ERHO")
			return Channel::ERHO;
		else if (chan == "MRHO")
			return Channel::MRHO;
		else if (chan == "EKAx")
			return Channel::EKAx;
		else if (chan == "MKAx")
			return Channel::MKAx;
		else if (chan == "nOMEGA")
			return Channel::nOMEGA;
		else if (chan == "nETA")
			return Channel::nETA;
		else if (chan == "nETAi")
			return Channel::nETAi;
		else if (chan == "nPHI")
			return Channel::nPHI;
		else if (chan == "EDs")
			return Channel::EDs;
		//else if (chan == "ExpALL")
			//return Channel::ExpALL;
		else
			throw std::invalid_argument("Decay channel "
						+ chan + " is unknwon");
	}

	std::string toString(Channel chan) {
		switch (chan) {
			case Channel::nnn:
				return "nnn";
			case Channel::nGAMMA:
				return "nGAMMA";
			case Channel::nEE:
				return "nEE";
			case Channel::nEE_o:
				return "nEE_o";
			case Channel::nEM:
				return "nEM";
			case Channel::nEM_o:
				return "nEM_o";
			case Channel::nMM:
				return "nMM";
			case Channel::nMM_o:
				return "nMM_o";
			case Channel::nET:
				return "nET";
			case Channel::nET_o:
				return "nET_o";
			case Channel::nMT:
				return "nMT";
			case Channel::nMT_o:
				return "nMT_o";
			case Channel::nPI0:
				return "nPI0";
			case Channel::EPI:
				return "EPI";
			case Channel::MPI:
				return "MPI";
			case Channel::TPI:
				return "TPI";
			case Channel::EKA:
				return "EKA";
			case Channel::MKA:
				return "MKA";
			case Channel::nRHO0:
				return "nRHO0";
			case Channel::ERHO:
				return "ERHO";
			case Channel::MRHO:
				return "MRHO";
			case Channel::EKAx:
				return "EKAx";
			case Channel::MKAx:
				return "MKAx";
			case Channel::nOMEGA:
				return "nOMEGA";
			case Channel::nETA:
				return "nETA";
			case Channel::nETAi:
				return "nETAi";
			case Channel::nPHI:
				return "nPHI";
			case Channel::EDs:
				return "EDs";
			default:
				throw std::invalid_argument("Decay channel "
							+ toString(chan) + " is unknwon");
		}
	}

	std::vector<Channel> Channels() {
		return { Channel::nnn,
			 Channel::nGAMMA,
			 Channel::nEE,
			 Channel::nEM,
			 Channel::nMM,
			 Channel::nET,
			 Channel::nMT,
			 Channel::nPI0,
			 Channel::EPI,
			 Channel::MPI,
			 Channel::TPI,
			 Channel::EKA,
			 Channel::MKA,
			 Channel::nRHO0,
			 Channel::ERHO,
			 Channel::MRHO,
			 Channel::EKAx,
			 Channel::MKAx,
			 Channel::nOMEGA,	
			 Channel::nETA,
			 Channel::nETAi,
			 Channel::nPHI,
			 Channel::EDs };
	}

	// these correspons to exp all
	std::vector<Channel> Detections() {
		return { Channel::nEE,
			 Channel::nEM,
			 Channel::nMM,
			 Channel::nPI0,
			 Channel::EPI,
			 Channel::MPI };
	}


	std::vector<double> Masses(Channel chan) {
		switch(chan) {
			//DECAYS
			//case Channel::ALL:
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
			case Channel::EDs:
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
				return {Const::MMuon, Const::MKaonx};
			default:
				throw std::invalid_argument("Decay channel "
							+ toString(chan) + " is unknwon");
		}
	}

	std::vector<int> Pdgs(Channel chan) {
		switch(chan) {
			//DECAYS
			//case Channel::ALL:
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
			case Channel::EDs:
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
			default:
				throw std::invalid_argument("Decay channel "
							+ toString(chan) + " is unknwon");
		}
	}
}
