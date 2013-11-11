#pragma once
#include "stdafx.h"
#include "armadillo"

using namespace arma;

class Scheme
{
public:
	Scheme(void);
	~Scheme(void);
	mat Scheme::explicitScheme(int nSteps, double tSteps, double(*u_s)(double));
	mat Scheme::implicitScheme(int nSteps, double tSteps, double(*u_s)(double));
	mat Scheme::crankNicolsonScheme(int nSteps, double tSteps, double(*u_s)(double));
};

