#include <iostream>
#include <string>

int main()
{
	std::string pt = "this/is/a/path/with.extension";
	std::string ot = "this/is/a/path/without";
	std::cout << pt << std::endl;
	std::cout << ot << std::endl;
	if (pt.find_last_of('.') != std::string::npos)
		pt.insert(pt.find_last_of('.'), "_a_nice");
	else
		pt.append("_a_nice");

	if (ot.find_last_of('.') != std::string::npos)
		ot.insert(ot.find_last_of('.'), "_a_nice_extension");
	else
		ot.append("_a_nice_extension");
	std::cout << std::endl;
	std::cout << pt << std::endl;
	std::cout << ot << std::endl;
	return 0;
}
