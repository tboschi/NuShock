/*
 * Detector class
 * Collections of relevant properties of near detector
 *
 * Author: Tommaso Boschi
 */

#ifndef TRACKER_H
#define TRACKER_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>

#include "tools/RNG.h"

#include "physics/Const.h"
#include "physics/Particle.h"

#include "detector/Detector.h"
#include "detector/Track.h"

class Tracker : public Detector
{
	public:
		using Event = std::pair<Particle, Track>;

	public:
		Tracker(std::string det, std::string track);
		Tracker(std::string card);

		Event GenerateEvent(Particle &&P) const;
		Event GenerateEvent(std::string mod, Particle &&P) const;
		bool Reconstruct(Event &evt) const;
		void ComputeTrack(Event &evt) const;
		//double ComputeSagitta(Event &evt) const;
		bool DetectSagitta(Event &evt) const;
		void Smearing(Event &evt) const;
		bool IsDecayed(const Particle &P, double dx) const;
		bool IsDetectable(const Event &evt) const;
		bool IsDetectable(const Track &T) const;
		bool IsDetectable(const Particle &P) const;

		std::array<Event, 2> Pi0Decay(Event &&pi0) const;
		std::vector<Event> MisIdentify(std::vector<Tracker::Event> &&evts) const;

	private:
		void Init(const CardDealer &cd);

		std::map<std::string, double> _thresholds;
		std::map<std::string, double> _resolutions;
		std::map<std::string, double> _identifies;
		double _step;
		bool kVerbose;
};

#endif
