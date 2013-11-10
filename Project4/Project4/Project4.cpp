// Project4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;
using namespace arma;

double u(double x);

int _tmain(int argc, _TCHAR* argv[])
{
	int iTimeSteps = 100; 
	int iNbSteps;
	bool bWantToQuit = false;

	while (!bWantToQuit)
	{
		cout << "Enter the number of steps you want \t";
		cin >> iNbSteps;
		fflush(stdin);
		int nScheme = 0;
		Scheme scheme = Scheme();
		cout << "To choose your algo, press 1 for the Explicit Scheme \n 2 for the Implicit Scheme \n and 3 for the Crank Nicholson Scheme \n";
		cout << "Any other key will quit \n";
		cin >> nScheme;
		switch (nScheme)
		{
		case 1:
			scheme.explicitScheme(iNbSteps, iTimeSteps, u);
			break;
		case 2:
			scheme.implicitScheme(iNbSteps, iTimeSteps, u);
			break;
		case 3:
			scheme.crankNicholsonScheme(iNbSteps, iTimeSteps, u);
			break;
		default:
			printf("You shall not be there  !! \n");
			//printf("(To choose your algo, press 1 for the Explicit Scheme \n 2 for the Implicit Scheme \n and 3 for the Crank Nicholson Scheme \n");
			bWantToQuit = true;
			break;
		}
	}

	return 0;
}

double u(double x)
{
	return 1- x;
}

