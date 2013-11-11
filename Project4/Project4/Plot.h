#pragma once
#include "armadillo"
#include "stdafx.h"

using namespace arma;
class Plot
{
public:
	Plot();
	~Plot();
	void WriteToFile(vec M, string path);
};