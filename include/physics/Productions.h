#ifndef PRODUCTIONS_H
#define PRODUCTIONS_H

#include <stdexcept>
#include <string>
#include <vector>

#include "physics/Const.h"

namespace Production
{
	enum class Channel {
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

	Channel fromString(std::string chan);
	std::string toString(Channel chan);

	// return all channels
	std::vector<Channel> Channels();
	std::vector<double> Masses(Channel chan);
	std::vector<int> Pdgs(Channel chan);
	unsigned int Ns(Channel chan);
}

inline std::ostream & operator<<(std::ostream &os, const Production::Channel &chan) {
	return os << Production::toString(chan);
}

#endif
