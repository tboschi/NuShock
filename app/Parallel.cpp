#include <iostream>

#include "TVector3.h"
#include "TRotation.h"

int main()
{
	TVector3 z(-1, 5, 1);
	TVector3 u(1, 2, 3);
	TRotation a;
	a.SetZAxis(u);

	z.Print();
	z = a * z;
	z.Print();
	return 0;
}
