#pragma once
#include "stdafx.h"
#include "armadillo"

using namespace arma;

class Scheme
{
public:
	Scheme(void);
	~Scheme(void);
	vec Scheme::explicitScheme(int nSteps, double tSteps, double(*u_s)(double));
	vec Scheme::implicitScheme(int nSteps, double tSteps, double(*u_s)(double));
	vec Scheme::crankNicholsonScheme(int nSteps, double tSteps, double(*u_s)(double));
};

