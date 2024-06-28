#ifndef OPTIMIZE_R
#define OPTIMIZE_R

#include <iostream>
#include <cmath>
#include <cfloat>

// Optim class: virtual class to find max/min of univariate function in R manner
// Usually a subclass is used to define a substantiated function
// Member function (public): value, evaluate function with doulbe x
// Other parameters may be added into the subclass as members
class Optim
{
public:
	virtual double value(double x) = 0;//x = alpha
	virtual ~Optim() {}
};

// Define optimize function
double optimize(Optim* optim, double lower, double upper, bool maximum, double tol);
#endif