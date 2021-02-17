#include<iostream>
#include <iomanip>
#include "Math/Numerical.hpp"

using namespace Panda;
using namespace std;

int main(int, char**)
{
    cout.precision(15);
    cout.setf(ios::fixed);

	NewtonRaphson<double, double>::nr_f f = [](double x) { return pow(x, 2.0) - 2*x + 1; };
	NewtonRaphson<double, double>::nr_fprime fprime = [](double x) {return 2 * pow(x, 1.0) - 2.0; };
    auto result = NewtonRaphson<double, double>::Solve(2, f, fprime);
    cout << "root of equation x^2 - 2x + 1 = 0 is: " << result << endl;

	getchar();

    return 0;
}