#ifndef CHANNELS_H
#define CHANNELS_H

#include <string>
#include <vector>

struct Channel
{
	enum class Type
	{
		no_type = -1,
		decayrates,
		production 
	};

	enum Name {
		undefined = -1,

		ALL,

		//decay modes
		//unclassified
		nnn,		//3 body	N -> 3 nu
		nGAMMA,	//2 body	N -> nu photon
		//pure leptonic
		nEE,		//3 body	N -> nu e e
		nEE_o,
		nEM,		//3 body	N -> nu e mu (via U_e)
		nEM_o,
		nMM,		//3 body	N -> nu mu mu
	        nMM_o,	
		nET,		//3 body	N -> nu e tau (via U_e)
	        nET_o,	
		nMT,		//3 body	N -> nu mu tau (via U_m)
	        nMT_o,	
		//pion
		nPI0,		//2 body	N -> nu pi0
		EPI,		//2 body	N -> e pi
		MPI,		//2 body	N -> mu pi
		TPI,		//2 body	N -> tau pi
		//kaon
		EKA,		//2 body	N -> e K
		MKA,		//2 body	N -> mu K
		//rho decay 100% in pions
		nRHO0,		//2 body	N -> rho0
		ERHO,		//2 body	N -> e rho
		MRHO,		//2 body	M -> mu rho
		//kaon*
		EKAx,		//2 body	N -> e K*
		MKAx,		//2 body	N -> mu K*
		//other (eta, phi, omega.. )
		nOMEGA,	//2 body	N -> nu w
		nETA,		//2 body	N -> nu eta
		nETAi,		//2 body	N -> nu eta'
		nPHI,		//2 body	N -> nu phi
		//charm
		ECHARM,	//2 body	N -> e D+
		//Channels for experimental comparison (EPI, MPI, nEE, nEM, nMM)
		ExpALL,	//

		//production modes
		//pure leptonic
		MuonE,		//3 body	mu  -> nu N e	(via Ue)
		MuonM,		//3 body	mu  -> nu N e	(via Um)
		TauEE,		//3 body	tau -> nu N e	(via Ue)
		TauET,		//3 body	tau -> nu N e	(via Ut)
		TauMM,		//3 body	tau -> nu N mu	(via Um)
		TauMT,		//3 body	tau -> nu N mu	(via Ut)
		TauPI,		//2 body	tau -> N pi	(via Ut)
		Tau2PI,	//3 body	tau -> N pi pi	(via Ut)
		//pseudomeson leptonic
		PionE,		//2 body	pi -> N e	(via Ue)
		PionM,		//2 body	pi -> N mu	(via Um)
		KaonE,		//2 body	K  -> N e	(via Ue)
		KaonM,		//2 body	K  -> N mu	(via Um)
		CharmE,	//2 body	Ds -> N e	(via Ue)
		CharmM,	//2 body	Ds -> N mu	(via Um)
		CharmT,	//2 body	Ds -> N tau	(via Ut)
		//pseudomeson semileptonic
		Kaon0E,	//3 body	K0 -> pi+ N e	(via Ue)
		Kaon0M,	//3 body	K0 -> pi+ N mu	(via Um)
		KaonCE,	//3 body	K+ -> pi0 N e	(via Ue)
		KaonCM 	//3 body	K+ -> pi0 N mu	(via Um)
	};

	static Name fromString(const std::string &chan) {
		if (chan == "ALL")			//1
			return Name::ALL;
		//decay channels
		else if (chan == "nnn")			//2
			return Name::nnn;
		else if (chan == "nGAMMA")
			return Name::nGAMMA;
		else if (chan == "nEE")
			return Name::nEE;
		else if (chan == "nEM")
			return Name::nEM;
		else if (chan == "nMM")
			return Name::nMM;
		else if (chan == "nET")
			return Name::nET;
		else if (chan == "nMT")
			return Name::nMT;
		else if (chan == "nPI0")
			return Name::nPI0;
		else if (chan == "EPI")
			return Name::EPI;
		else if (chan == "MPI")
			return Name::MPI;
		else if (chan == "TPI")
			return Name::TPI;
		else if (chan == "EKA")
			return Name::EKA;
		else if (chan == "MKA")
			return Name::MKA;
		else if (chan == "nRHO0")
			return Name::nRHO0;
		else if (chan == "ERHO")
			return Name::ERHO;
		else if (chan == "MRHO")
			return Name::MRHO;
		else if (chan == "EKAx")
			return Name::EKAx;
		else if (chan == "MKAx")
			return Name::MKAx;
		else if (chan == "nOMEGA")
			return Name::nOMEGA;
		else if (chan == "nETA")
			return Name::nETA;
		else if (chan == "nETAi")
			return Name::nETAi;
		else if (chan == "nPHI")
			return Name::nPHI;
		else if (chan == "ECHARM")
			return Name::ECHARM;
		else if (chan == "ExpALL")		//25
			return Name::ExpALL;
		//production channels
		else if (chan == "MuonE")		//26
			return Name::MuonE;
		else if (chan == "MuonM")
			return Name::MuonM;
		else if (chan == "TauEE")
			return Name::TauEE;
		else if (chan == "TauET")
			return Name::TauET;
		else if (chan == "TauMM")
			return Name::TauMM;
		else if (chan == "TauMT")
			return Name::TauMT;
		else if (chan == "TauPI")
			return Name::TauPI;
		else if (chan == "Tau2PI")
			return Name::Tau2PI;
		else if (chan == "PionE")
			return Name::PionE;
		else if (chan == "PionM")
			return Name::PionM;
		else if (chan == "KaonE")
			return Name::KaonE;
		else if (chan == "KaonM")
			return Name::KaonM;
		else if (chan == "CharmE")
			return Name::CharmE;
		else if (chan == "CharmM")
			return Name::CharmM;
		else if (chan == "CharmT")
			return Name::CharmT;
		else if (chan == "Kaon0E")
			return Name::Kaon0E;
		else if (chan == "Kaon0M")
			return Name::Kaon0M;
		else if (chan == "KaonCE")
			return Name::KaonCE;
		else if (chan == "KaonCM")
			return Name::KaonCM;
		else
			return Name::undefined;
	}

	static std::string toString(const Name &chan) {
		switch (chan) {
			case Name::undefined:
				return "undefined";
			case Name::ALL:
				return "ALL";	
			case Name::nnn:
				return "nnn";
			case Name::nGAMMA:
				return "nGAMMA";
			case Name::nEE:
				return "nEE";
			case Name::nEE_o:
				return "nEE_o";
			case Name::nEM:
				return "nEM";
			case Name::nEM_o:
				return "nEM_o";
			case Name::nMM:
				return "nMM";
			case Name::nMM_o:
				return "nMM_o";
			case Name::nET:
				return "nET";
			case Name::nET_o:
				return "nET_o";
			case Name::nMT:
				return "nMT";
			case Name::nMT_o:
				return "nMT_o";
			case Name::nPI0:
				return "nPI0";
			case Name::EPI:
				return "EPI";
			case Name::MPI:
				return "MPI";
			case Name::TPI:
				return "TPI";
			case Name::EKA:
				return "EKA";
			case Name::MKA:
				return "MKA";
			case Name::nRHO0:
				return "nRHO0";
			case Name::ERHO:
				return "ERHO";
			case Name::MRHO:
				return "MRHO";
			case Name::EKAx:
				return "EKAx";
			case Name::MKAx:
				return "MKAx";
			case Name::nOMEGA:
				return "nOMEGA";
			case Name::nETA:
				return "nETA";
			case Name::nETAi:
				return "nETAi";
			case Name::nPHI:
				return "nPHI";
			case Name::ECHARM:
				return "ECHARM";
			case Name::ExpALL:
				return "ExpALL";		//25
			case Name::MuonE:
				return "MuonE";		//26
			case Name::MuonM:
				return "MuonM";
			case Name::TauEE:
				return "TauEE";
			case Name::TauET:
				return "TauET";
			case Name::TauMM:
				return "TauMM";
			case Name::TauMT:
				return "TauMT";
			case Name::TauPI:
				return "TauPI";
			case Name::Tau2PI:
				return "Tau2PI";
			case Name::PionE:
				return "PionE";
			case Name::PionM:
				return "PionM";
			case Name::KaonE:
				return "KaonE";
			case Name::KaonM:
				return "KaonM";
			case Name::CharmE:
				return "CharmE";
			case Name::CharmM:
				return "CharmM";
			case Name::CharmT:
				return "CharmT";
			case Name::Kaon0E:
				return "Kaon0E";
			case Name::Kaon0M:
				return "Kaon0M";
			case Name::KaonCE:
				return "KaonCE";
			case Name::KaonCM:
				return "KaonCM";
			default:
				return "undefined";
		}
	}

	static Type whichType(const Name &chan) {
		if (chan <= Name::ALL)
			return Type::no_type;
		if (chan <= Name::ExpALL)
			return Type::decayrates;
		if (chan <= Name::KaonCM)
			return Type::production;
		return Type::no_type;
	}

	static std::vector<Channel::Name> Decays() {
		return { Channel::nnn,		//3 body	N -> 3 nu
			 Channel::nGAMMA,	//2 body	N -> nu photon
			 Channel::nEE,		//3 body	N -> nu e e
			 Channel::nEM,		//3 body	N -> nu e mu (via U_e)
			 Channel::nMM,		//3 body	N -> nu mu mu
			 Channel::nET,		//3 body	N -> nu e tau (via U_e)
			 Channel::nMT,		//3 body	N -> nu mu tau (via U_m)
			 Channel::nPI0,		//2 body	N -> nu pi0
			 Channel::EPI,		//2 body	N -> e pi
			 Channel::MPI,		//2 body	N -> mu pi
			 Channel::TPI,		//2 body	N -> tau pi
			 Channel::EKA,		//2 body	N -> e K
			 Channel::MKA,		//2 body	N -> mu K
			 Channel::nRHO0,		//2 body	N -> rho0
			 Channel::ERHO,		//2 body	N -> e rho
			 Channel::MRHO,		//2 body	M -> mu rho
			 Channel::EKAx,		//2 body	N -> e K*
			 Channel::MKAx,		//2 body	N -> mu K*
			 Channel::nOMEGA,	//2 body	N -> nu w
			 Channel::nETA,		//2 body	N -> nu eta
			 Channel::nETAi,		//2 body	N -> nu eta'
			 Channel::nPHI,		//2 body	N -> nu phi
			 Channel::ECHARM };
	}

	static std::vector<Channel::Name> Detections() {
		return { Channel::nEE,		//3 body	N -> nu e e
			 Channel::nEM,		//3 body	N -> nu e mu (via U_e)
			 Channel::nMM,		//3 body	N -> nu mu mu
			 Channel::nPI0,		//2 body	N -> nu pi0
			 Channel::EPI,		//2 body	N -> e pi
			 Channel::MPI };		//2 body	N -> mu pi
	}

	static std::vector<Channel::Name> Productions() {
		return { Channel::MuonE,	//3 body	mu  -> nu N e	(via Ue)
			 Channel::MuonM,	//3 body	mu  -> nu N e	(via Um)
			 Channel::TauEE,	//3 body	tau -> nu N e	(via Ue)
			 Channel::TauET,	//3 body	tau -> nu N e	(via Ut)
			 Channel::TauMM,	//3 body	tau -> nu N mu	(via Um)
			 Channel::TauMT,	//3 body	tau -> nu N mu	(via Ut)
			 Channel::TauPI,	//2 body	tau -> N pi	(via Ut)
			 Channel::Tau2PI,	//3 body	tau -> N pi pi	(via Ut)
			 Channel::PionE,	//2 body	pi -> N e	(via Ue)
			 Channel::PionM,	//2 body	pi -> N mu	(via Um)
			 Channel::KaonE,	//2 body	K  -> N e	(via Ue)
			 Channel::KaonM,	//2 body	K  -> N mu	(via Um)
			 Channel::CharmE,	//2 body	Ds -> N e	(via Ue)
			 Channel::CharmM,	//2 body	Ds -> N mu	(via Um)
			 Channel::CharmT,	//2 body	Ds -> N tau	(via Ut)
			 Channel::Kaon0E,	//3 body	K0 -> pi+ N e	(via Ue)
			 Channel::Kaon0M,	//3 body	K0 -> pi+ N mu	(via Um)
			 Channel::KaonCE,	//3 body	K+ -> pi0 N e	(via Ue)
			 Channel::KaonCM };
	}
};

inline std::ostream & operator<<(std::ostream &os, const Channel::Name &chan) {
	return os << Channel::toString(chan);
}
#endif
