#include "Detector.h"

Detector::Detector(std::string ConfigName)
{
	std::fopen ConfigFile(ConfigName.c_str());
	std::string Line, Key;
	double Element;
	EnergyEfficiency EnEff;
	while (getline(ConfigFile, Line))
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
			mapDetector[Key].push_ = Element;
		}
	 }
}

void ListKey()
{
	std::map<std::string, double>::iterator it;
	for (it = mapDetector.begin(); it != mapDetector.end(); ++it)
		std::cout << it->frist << std::endl;
}

double GetElement(std::string Key)
{
	return mapDetector[Key];
}

double Efficiency(std::string Channel, double Energy)
{
	double diff = Energy;
	std::vector<EnergyEfficiency>::iterator it, ip;
	for (it = mapEfficiency[Channel].begin(); it != mapEfficiency[Channel].end(); ++it)
	{
		if (diff > (Energy - it->second.E))
		{
			diff = Energy - it->second.E;
			ip = it;
		}
	}

	return Energy*ip->second.f/sqrt(ip->second.E);
}
