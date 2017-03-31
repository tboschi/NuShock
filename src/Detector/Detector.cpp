#include "Detector.h"

Detector::Detector(std::string ConfigName)
{
	std::ifstream ConfigFile(ConfigName.c_str());

	std::string Line, Key;
	std::stringstream ssL;
	double Element;
	EnergyEfficiency EnEff;

	while (std::getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		if (Line[0] == 'W')
		{
			ssL >> Key >> EnEff.E >> EnEff.f;
			Key.erase(Key.begin());
			mapEfficiency[Key].push_back(EnEff);
		}
		else
		{
			ssL >> Key >> Element;
			mapDetector[Key] = Element;
		}
	 }
}

std::vector<std::string> Detector::ListKey()
{
	std::vector<std::string> List;
	std::map<std::string, double>::iterator it;
	for (it = mapDetector.begin(); it != mapDetector.end(); ++it)
		List.push_back(it->first);
	return List;
}

std::vector<std::string> Detector::ListChannel()
{
	std::vector<std::string> List;
	std::map<std::string, std::vector<EnergyEfficiency> >::iterator it;
	for (it = mapEfficiency.begin(); it != mapEfficiency.end(); ++it)
		List.push_back(it->first);
	return List;
}

double Detector::GetElement(std::string Key)
{
	return mapDetector[Key];
}

double Detector::Efficiency(std::string Channel, double Energy)
{
	double diff = 10.0, sdiff = 10.0;
	double EA = 0, EB = 0;
	double fA, fB;


	std::vector<EnergyEfficiency>::iterator it, iA, iB;
	for (it = mapEfficiency[Channel].begin(); it != mapEfficiency[Channel].end(); ++it)
	{
		if (Energy == it->E)
			return it->f;
		else if (Energy > it->E)
		{
			if (diff > (Energy - it->E))
			{
				diff = Energy - it->E;
				EA = it->E;
				fA = it->f;
			}
		}
		else 
		{
			if (sdiff > (it->E - Energy))
			{
				sdiff = it->E - Energy;
				EB = it->E;
				fB = it->f;
			}
		}
	}

	if (EA == 0)
		return fB;
	else if (EB == 0)
		return fA;
	else
	{
		double mFactor = (fB - fA)/(EB - EA);
		double qFactor = fA - mFactor * EA;
	
		return Energy*mFactor+qFactor;
	}
}

double Detector::EnergySigma(std::string Channel, double Energy)
{
	return Energy*(Energy*Efficiency(Channel, Energy))/sqrt(Energy);
}
