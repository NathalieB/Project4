// Project4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;
using namespace arma;

double u(double x);
void WriteToFile(mat M, string path);

int _tmain(int argc, _TCHAR* argv[])
{
	int iNbSteps = 5;
	int truc = pow(iNbSteps + 1, 2);
	int iTimeSteps = (float)truc/0.4;
	bool bWantToQuit = false;

	while(!bWantToQuit)
	{
		printf("The configuration is the following : \n \t time steps : %d \t spatial steps: %d \n", iTimeSteps,iNbSteps);
		int nScheme = 0;
		Scheme scheme = Scheme();
		cout << "To choose your algo, press 1 for the Explicit Scheme \n 2 for the Implicit Scheme \n and 3 for the Crank Nicolson Scheme \n";
		cout << "Any other key will quit \n";
		cin >> nScheme;
		switch(nScheme)
		{
		case 1:
			WriteToFile(scheme.explicitScheme(iNbSteps, iTimeSteps, u),"test.txt");
			break;
		case 2:
			WriteToFile(scheme.implicitScheme(iNbSteps, iTimeSteps, u), "test2.txt");
			break;
		case 3:
			WriteToFile(scheme.crankNicolsonScheme(iNbSteps, iTimeSteps, u),"testCN.txt");
			break;
		default:
			printf("You shall not be there  !! \n");
			//printf("(To choose your algo, press 1 for the Explicit Scheme \n 2 for the Implicit Scheme \n and 3 for the Crank nicolson Scheme \n");
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

// Write the content of a matrix to file for plotting in Matlab
// (now redundant because armadillo does the same thing)
void WriteToFile(mat M, string path)
{
	int nRow = M.n_rows;
	int nCol = M.n_cols;
	ofstream outfile; // file stream
	outfile.open(path); // attempt to open file for writing

	for (int i = 0; i < nRow; i++) // rows
	{
		for (int j = 0; j < nCol; j++) // columns
		{
			if (j > 0) // no leading space
			{
				outfile << ' ';
			}

			outfile << M(i,j);
		}

		outfile << endl; // newline character
	}

	outfile.close();
}