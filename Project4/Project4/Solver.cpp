#include "stdafx.h"
#include "Solver.h"

/*
	This class is a generic class for the:
		* Jacobi solver
	algorithms
*/

Solver::Solver()
{
}


Solver::~Solver()
{
}

// Function used to do the first step of the Gaussion elimination: the Forward Substitution
// We'll kill every term which prevents the matrix A to be an upper triangular one.
// Note the &'s needed to make Armadillo objects be passed by reference!
bool Solver::forwardSubstitutionVector(vec &v_a, vec &v_d, vec &v_f)
{
	int sizeVector = v_f.n_rows;
	double coef = 0.0;
	// And then, we compute the forward substitution:
	for (int i= 1; i< sizeVector; i++)
	{
		coef = v_a[i] / v_d[i - 1];
		v_d[i] -= coef * v_a[i];
		v_a[i] -= coef*v_d[i-1];  // This operation just allows us to check if  every term of the first diago is null after the substitution
		v_f[i] -= coef*v_f[i-1];
		printf("a%d : %f \t b%d : %f \n",i, v_a[i],i, v_d[i]); // debug
	}

	//for (int i=0; i< sizeVector; i++) // debug
		//printf(" After \tf%d : %f \n", i, v_f[i]); // debug

	return true;
}

// Perform backward substitution
// Note the &'s needed to make Armadillo objects be passed by reference!
bool Solver::backwardSubstitutionVector(vec &v_f, vec &v_b, vec v_c, vec &v_Solution)
{
	int sizeVector = v_f.n_rows;

	// We save the last term:
	v_Solution[sizeVector-1] = v_f[sizeVector-1]/v_b[sizeVector-1];
	// and then we compute what's left
	for (int i = sizeVector-1; i > 0; i--)
	{
		v_Solution[i-1] = (v_f[i-1] - v_c[i-1]*v_Solution[i]) / v_b[i-1];
		//printf(" u%d : %f | ", i-1, v_Solution[i-1]); // Not displaying it... Pretty clear for little number of Row, but not in other cases
	}
	//printf("\n"); 
	return true;
}

// a: diag. elem | b: non diag. elem
void Solver::tridiagonalSolver(double a,double b,vec& v,vec& vNext)
{

	int sizeVector = v.n_rows;
	// First, we need to duplicate v_a, since we have to "nearly diagonals"
	vec v_a(sizeVector), v_c(sizeVector);
	vec v_d(sizeVector);
	// Initialization of our 2nd vector
	for (int i = 0; i < sizeVector; i++)
		v_a[i] = v_c[i] = b;

	for (int i = 0; i< sizeVector; i++)
		v_d[i] = a;

	// Calling Forward then Backward 
	forwardSubstitutionVector(v_a, v_d,vNext);
	backwardSubstitutionVector(vNext, v_d, v_c, v);
	
	return;
}