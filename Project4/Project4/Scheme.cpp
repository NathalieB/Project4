#include "stdafx.h"
#include "Scheme.h"
#include "Solver.h"
#include "Plot.h"

/*
	This class contains:
		* the Explicit Scheme
		* the Implicit Scheme
		* the Crank nicolson Scheme
	algorithms
*/

Scheme::Scheme(void)
{
}


Scheme::~Scheme(void)
{
}

mat Scheme::explicitScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps , time, u_s(x)

	double deltaX = 1.0 / (float)(nSteps + 1);
	double alpha = 0.4;
	double deltaT = alpha * pow(deltaX, 2);
	//double 
	tSteps = 1.0 / deltaT;

	vec v(nSteps), vNext(nSteps), u(nSteps);
	int printIndex = tSteps / 10;
	mat M(printIndex, nSteps);
	int index = 0;

	//Initialization 
	for (int i = 1; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		v(i) = -u_s(x);
	}

	// Boundary conditions
	v[0] = vNext[0] = 0;
	v[nSteps] = vNext[nSteps] = 0;
	//
	for(int t = 1; t <= tSteps; t ++)
	{
		if ((t - 1) % 204 == 0)
		{
			for (int i = 0; i < nSteps; i++)
				printf("%f \t", v(i));
			printf("\n");
		}

		for(int i = 1; i < nSteps; i++ )
		{
			vNext[i] = (1 - 2 * alpha) * v[i];
			double truc = vNext[i]; // XX
			if( i > 0 ) 
				vNext[i] += alpha * v[i-1]; // else += 0
			if( i < nSteps ) 
				vNext[i] += alpha * v[i+1]; // else += 0
		}
		v = vNext;

		/*if ((t - 1) % 204 == 0)
		{
			for (int i = 0; i < nSteps; i++)
				printf("%f \t", v(i));
			printf("\n");
		}*/

		for (int i = 0; i < nSteps; i++)
		{
			double x = (i + 1) * deltaX;
			u[i] = v[i] + u_s(x);
		}


		// Print to file !
		if (t % 10 == 0)
		{
			for (int i = 0; i < nSteps; i++)
				M(index, i) = u[i];
			index++;
		}
	}
	
	for (int i = 0; i < printIndex; i++)
	{
		for (int j = 0; j < nSteps; j++)
			printf("%f \t", M(i, j));
	}

	return M;
}

mat Scheme::implicitScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps (# of interior points), time, tSteps, u_s(x)

	double deltaX = 1.0 / (nSteps + 1);
	double deltaT = 1.0 / tSteps;
	double alpha = deltaT / pow(deltaX,2);

	vec v(nSteps), vNext(nSteps), u(nSteps);
	int printIndex = tSteps / 10;
	mat M(printIndex, nSteps);
	int index = 0;

	for( int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		vNext[i] = v[i] = -u_s(x);
	}
	// Boundary conditions
	v[0] = vNext[0] = v[nSteps] = vNext[nSteps] = 0;
	
	double a = -alpha;        // diagonal element
	double b = 1 + 2*alpha;   // off-diagonal element

	Solver solv = Solver();
	for(int t = 1; t <= tSteps; t++)
	{
		solv.tridiagonalSolver(a, b, v, vNext);
		
		// Initial conditions
		vNext[0] = vNext[nSteps] = 0;
		v = vNext;

		for (int i = 0; i < nSteps; i++)
		{
			double x = (i + 1) * deltaX;
			u[i] = v[i] + u_s(x);
		}

		// Print to file !
		if (t % 10 == 0)
		{
			for (int i = 0; i < nSteps; i++)
				M(index, i) = u[i];
			index++;
		}
	}

	for (int i = 0; i < printIndex; i++)
	{
		for (int j = 0; j < nSteps; j++)
			printf("%f \t", M(i, j));
	}

	return M;	
}

mat Scheme::crankNicolsonScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps (# of interior points), time, tSteps, u_s(x)

	double deltaX = 1.0 / (nSteps + 1);
	double deltaT = 1.0 / tSteps;
	double alpha = deltaT / pow(deltaX,2);

	vec v(nSteps), w(nSteps), u(nSteps);

	int printIndex = tSteps / 10;
	mat M(printIndex, nSteps);
	int index = 0;

	for(int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		v[i] = -u_s(x);
	}


	// Boundary conditions
	v[0] = w[0] = v[nSteps] = w[nSteps] = 0;


	double a = 2 * (1 + alpha);    // diagonal element
	double b = -alpha;             // off-diagonal element
	Solver solv = Solver();
	for(int t = 1; t <= tSteps; t++)
	{
		for( int i = 0; i < nSteps; i++ )
		{
			w[i] = 2 * (1 - alpha) * v[i];
			if( i > 0 ) 
				w[i] += alpha * v[i - 1]; // else += 0
			if( i < nSteps - 1 ) 
				w[i] += alpha * v[i + 1]; // else += 0
		}	
		solv.tridiagonalSolver(a, b, w, v);

		
		// Boundary conditions:
		v[0] = v[nSteps] = 0;

		for (int i = 0; i < nSteps; i++)
		{
			double x = (i + 1) * deltaX;
			u[i] = v[i] + u_s(x);
		}

		// Print to file !
		if (t % 10 == 0)
		{
			for (int i = 0; i < nSteps; i++)
				M(index, i) = u[i];
			index++;
		}
	}

	for (int i = 0; i < printIndex; i++)
	{
		for (int j = 0; j < nSteps; j++)
			printf("%f \t", M(i, j));
	}

	return M;
}