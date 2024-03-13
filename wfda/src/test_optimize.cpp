// Unit test
#include "optimize.h"
using namespace std;

// Define a class of hyperbolic function: f(x) = a * x^2 + b * x + c
class Parabol: public Optim
{
private:
	const double a;
	const double b;
	const double c;
public:
	Parabol(double a_, double b_, double c_) : a(a_), b(b_), c(c_) {}
	double value(double x) {
		return a * pow(x, 2) + b * x + c;
	}
};

// Main function to minimize a hyperbolic function
int main() {
	Parabol parabol1(1, -5, 3);
	double x_min = optimize(&parabol1, 0, 5, false, 1e-3);
	cout << x_min << endl;
	Parabol parabol2(-1, -5, 3);
	double x_max = optimize(&parabol2, -5, 0, true, 1e-3);
	cout << x_max << endl;
	return 0;
}