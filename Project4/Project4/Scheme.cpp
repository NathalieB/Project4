#include "stdafx.h"
#include "Scheme.h"
#include "Solver.h"

/*
	This class contains:
		* the Explicit Scheme
		* the Implicit Scheme
		* the Crank Nicholson Scheme
	algorithms
*/

Scheme::Scheme(void)
{
}


Scheme::~Scheme(void)
{
}

vec Scheme::explicitScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps , time, u_s(x)

	double deltaX = 1.0 / (float)(nSteps + 1);
	double alpha = 0.5;
	double deltaT = alpha * pow(deltaX,2);
	//double 
	tSteps = 1.0 / deltaT;

	vec v(nSteps), vNext(nSteps), u(nSteps); 

	for(int i= 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		v(i) = -u_s(x);
	}
		
	for(int t = 1; t <= tSteps; t ++)
	{
		for(int i = 1; i < nSteps; i++ )
		{
			vNext[i] = (1 - 2 * alpha) * v[i];
			if( i > 0 ) 
				vNext[i] += alpha * v[i-1]; // else += 0
			if( i < nSteps ) 
				vNext[i] += alpha * v[i+1]; // else += 0
		}
		v = vNext;
	}

	for(int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		u[i] = v[i] + u_s(x);
	}

	return u;
}

vec Scheme::implicitScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps (# of interior points), time, tSteps, u_s(x)

	double deltaX = 1.0 / (nSteps + 1);
	double deltaT = 1.0 / tSteps;
	double alpha = deltaT / pow(deltaX,2);

	vec v(nSteps), vNext(nSteps), u(nSteps);


	for( int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		v[i] = -u_s(x);
	}

	// Missing Boundary conditions ?
	double a = -alpha;        // diagonal element
	double b = 1 + 2*alpha;   // off-diagonal element

	Solver solv = Solver();
	for(int t = 1; t <= tSteps; t++)
	{
		solv.tridiagonalSolver(a, b, v, vNext);
		v = vNext; 
	}

	for(int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		u[i] = v[i] + u_s(x);
	}

	return u;	
}

vec Scheme::crankNicholsonScheme(int nSteps, double tSteps, double(*u_s)(double))
{
	//input: nSteps (# of interior points), time, tSteps, u_s(x)

	double deltaX = 1.0 / (nSteps + 1);
	double deltaT = 1.0 / tSteps;
	double alpha = deltaT / pow(deltaX,2);

	vec v(nSteps), w(nSteps), u(nSteps);

	for(int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		v[i] = -u_s(x);
	}

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
	}

	for(int i = 0; i < nSteps; i++)
	{
		double x = (i + 1) * deltaX;
		u[i] = v[i] + u_s(x);
	}

	return u;
}