#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		std::cout << "insert mu e bar components" << std::endl;
		return 1;
	}

	double m = strtod(argv[1], NULL);
	double e = strtod(argv[2], NULL);
	double b = strtod(argv[3], NULL);

	int b_m, b_e, b_b;
	std::cout << "m component: ";
	std::cin >> b_m;
	std::cout << "e component: ";
	std::cin >> b_e;
	std::cout << "b component: ";
	std::cin >> b_b;

	std::cout << std::endl;
	double sum = (b_m == 0 ? 0.3 : b_m) * m +
		     (b_e == 0 ? 0.3 : b_e) * e +
		     (b_b == 0 ? 0.3 : b_b) * b;

	std::cout << "Total " << int(sum) << std::endl;
	return 0;
}       
