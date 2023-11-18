#include <iostream>
#include <math.H>
using namespace std;

double function(double x);
double simpsonsIntegral(double a, double b, int n);
double trapezoidIntegral(double a, double b, int n);
double integralT(double a, double b, double eps, int n);
double integralS(double a, double b, double eps, int n);

void main() {
	double eps = 0.00001, a = 3, b = 4.254, integral2, integral1, R;
	int n = 2;
	int choice;
	cout << "Choose the way to find integral. N=1 - Trapezoid, N=2 - Simpson. N=";
	cin >> choice;
	cout << endl;
	if (choice == 1) {
		integral2 = integralT(a, b, eps, n);
		integral1 = trapezoidIntegral(a, b, n);
		R = 4 * (integral1 - integral2) / 3;
		cout << endl << "error = " << R;
	}
	if (choice == 2) {
		integral2 = integralS(a, b, eps, n);
		integral1 = simpsonsIntegral(a, b, n);
		R = (integral1 - integral2) * (pow(0.5, 4) - 1);
		cout << endl << "error = " << R;
	}
	if (choice != 2 && choice != 1) {
		cout << "wrong input";
	}
}

double function(double x) {
	return ((1 + pow(x, 2)) / (1 + pow(x, 3)));
}
double simpsonsIntegral(double a, double b, int n) {
	double h = (b - a) / n;
	double part1 = 0, part2 = 0;
	for (int i = 1; i < n; ++i) {
		if (i % 2 == 0)
			part1 += function(a + h * i);
		if (i % 2 == 1)
			part2 += function(a + (i + 1) * h);
	}
	return h / 3 * (function(a) + function(b) + 4 * part2 + 2 * part1);
}
double trapezoidIntegral(double a, double b, int n) {
	double h = (b - a) / n;
	double part = 0;
	for (int i = 1; i < n; ++i) {
		part += function(a + h * i);
	}

	return h / 2 * (function(a) + 2 * part + function(a + h * n));
}
double integralT(double a, double b, double eps, int n) {
	double Integr1, Integr2;
	Integr1 = trapezoidIntegral(a, b, n);
	do {
		Integr2 = Integr1;
		n = n * 2;
		Integr1 = trapezoidIntegral(a, b, n);
	} while (abs(Integr1 - Integr2) <= 3 * eps);
	cout << "IntegralT = " << Integr1;
	return Integr1;
}
double integralS(double a, double b, double eps, int n) {
	double Integr1, Integr2;
	Integr1 = simpsonsIntegral(a, b, n);
	do {
		Integr2 = Integr1;
		n = n * 2;
		Integr1 = simpsonsIntegral(a, b, n);
	} while (abs(Integr1 - Integr2) >= 15 * eps);
	cout << "Integral = " << Integr1;
	return Integr1;
}
