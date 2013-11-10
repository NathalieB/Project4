#pragma once
#include "stdafx.h"

#include <stdio.h>
#include <tchar.h>
#include "armadillo"

using namespace std;
using namespace arma;

class Solver
{
public:
	Solver();
	~Solver();
	bool Solver::forwardSubstitutionVector(vec &v_a, vec &v_d, vec &v_f);
	bool Solver::backwardSubstitutionVector(vec &v_f, vec&v_b, vec v_c, vec &v_Solution);
	void Solver::tridiagonalSolver(double a,double b,vec& v,vec& vNext); // a: diag. elem | b: non diag. elem
};

