#include "Detector.h"

Detector::Detector(std::string ConfigName)
{
	std::fopen ConfigFile(ConfigName.c_str());
	std::string Line, Key;
	double Element;
	while (getline(ConfigFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> Key >> Element;
		mapDetector[Key] = Element;
	}
}
