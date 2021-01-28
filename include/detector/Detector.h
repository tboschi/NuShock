/*
 * Detector class
 * Collections of relevant properties of near detector
 *
 * Author: Tommaso Boschi
 */

#ifndef DETECTOR_H
#define DETECTOR_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>

//ROOT include
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "tools/CardDealer.h"

#include "detector/Materials.h"
#include "detector/Track.h"
#include "physics/Const.h"
#include "physics/Channels.h"
#include "physics/Particle.h"

class Detector
{
	public:

		using Module = std::map<std::string, double>;

		Detector(const std::string &card);
		double Efficiency(Channel::Name chan, double mass, double energy) const;
		void LoadEfficiency(Channel::Name, std::string file);

		Material::Name MadeOf(std::string mod) const;

		double Exposure() const;
		double Exposure(std::string mod) const;
		double Scaling() const;
		double Scaling(std::string mod) const;
		double Weight() const;
		double Weight(std::string mod) const;
		double Baseline() const;
		double Baseline(std::string mod) const;
		double Length() const;
		double Length(std::string mod) const;
		double Height() const;
		double Height(std::string mod) const;
		double Width() const;
		double Width(std::string mod) const;
		double Section() const;
		double Section(std::string mod) const;
		double Radius() const;
		double Radius(std::string mod) const;
		double Volume() const;
		double Volume(std::string mod) const;
		double AngularAcceptance() const;
		double AngularAcceptance(std::string mod) const;
		bool AngularAccept(double theta, double phi) const;
		bool AngularAccept(std::string mod, double theta, double phi) const;
		bool AngularAccept(const Particle &P) const;
		bool AngularAccept(std::string mod, const Particle &P) const;

		double Probability(double tby) const;
		double Probability(std::string mod, double tby) const;

		std::vector<std::string> Modules() const;
		std::string WhichModule(const Track &t) const;
		bool IsInside(const Track &t) const;
		bool IsInside(std::string mod, const Track &t) const;
		bool IsContained(const Track &t) const;

		double MagneticField(std::string mod) const;

		double BeamEnergy() const;
		double POTs() const;
		double POTs(std::string hc) const;

		friend std::ostream& operator<<(std::ostream &os, const Detector &box);

	protected:
		std::unordered_map<std::string, Module> _modules;
		std::unordered_map<std::string, Material::Name> _materials;
		std::unordered_map<std::string, std::string> _shapes;
		std::unordered_map<std::string, double> _modes;
		//std::set<std::string> _efficiencies;
		std::string _default;

		// memoization
		mutable double _baseline, _length, _width, _height, _volume, _exposure, _weight;

		double _POTs, _Eb;

		
		//std::map<Channel::Name, std::shared_ptr<TH2D> > _mass_ener_func;
};

#endif
