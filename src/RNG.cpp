#include "tools/RNG.h"

namespace RNG {
	std::mt19937 _mt = std::mt19937(std::random_device()());
	//std::mt19937 _mt = std::mt19937(0);
}
