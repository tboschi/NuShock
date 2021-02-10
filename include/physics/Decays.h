#ifndef DECAYS_H
#define DECAYS_H

#include <string>
#include <vector>
#include <stdexcept>

#include "physics/Const.h"

namespace Decay
{
	enum class Channel {
		//ALL,

		//decay modes
		//unclassified
		nnn,		//3 body	N -> 3 nu
		nGAMMA,	//2 body	N -> nu photon
		//pure leptonic
		nEE,		//3 body	N -> nu e e
		nEE_o,		// special
		nEM,		//3 body	N -> nu e mu (via U_e)
		nEM_o,		// special
		nMM,		//3 body	N -> nu mu mu
	        nMM_o,		// special
		nET,		//3 body	N -> nu e tau (via U_e)
	        nET_o,		// special
		nMT,		//3 body	N -> nu mu tau (via U_m)
	        nMT_o,		// special
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
		EDs,	//2 body	N -> e D+
		//Channels for experimental comparison (EPI, MPI, nEE, nEM, nMM)
		ExpALL	//
	};


	Channel fromString(std::string chan);
	std::string toString(Channel chan);

	std::vector<Channel> Channels();
	std::vector<Channel> Detections();
	std::vector<double> Masses(Channel chan);
	std::vector<int> Pdgs(Channel chan);
}

inline std::ostream & operator<<(std::ostream &os, const Decay::Channel &chan) {
	return os << Decay::toString(chan);
}

#endif
