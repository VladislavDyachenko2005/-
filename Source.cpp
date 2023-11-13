#include <iostream>
#include <math.H>
using namespace std;

double function(double x);
double simpsonsIntegral(double a, double b,int n);

void main() {
	double eps = 0.00001, a=3, b=4.254, simps1, simps2;
	int n = 2;
	simps1 = simpsonsIntegral(a, b, n);
	do  {
		simps2 = simps1;
		n = n * 2;
		simps1= simpsonsIntegral(a, b, n);
	} while (abs(simps1 - simps2) >= 15*eps);
	cout << "Integral = " << simps1;
}

double function(double x) {
	return ((1 + pow(x, 2)) / (1 + pow(x, 3)));
}
double simpsonsIntegral(double a, double b, int n) {
	double h = (b-a) / n;
	double part1=0, part2=0;
	for (int i = 1; i <= n; i++) {
		if (n % 2 == 0) {
			part1 += function(a + h * i);
		}
		if (n % 2 == 1) {
			part2 += function(a + h * i);
		}
	}
	return h / 3 * (function(a) + function(b) + 4 * part1 + 2 * part2);
}