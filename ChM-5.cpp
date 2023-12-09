#include <iostream>
#include <math.H>
using namespace std;

double function(double x);
double function2(double x, double y);
double simpsonsIntegral(double a, double b, int n);
double trapezoidIntegral(double a, double b, int n);
double integralT(double a, double b, double eps, int n);
double integralS(double a, double b, double eps, int n);
double doubleIntegral(double nx, double ny, double lx, double ux, double ly, double uy);

void main() {
	double eps = 0.00001, a = 3, b = 4.254, integral2, integral1, R;
	int n = 2;
	int choice;
	cout << "Choose the way to find integral. N=1 - Trapezoid, N=2 - Simpson. N=";
	cin >> choice;
	cout << endl;
	if (choice == 1) {
		integral2 = integralT(a, b, eps, n);
		cout << "Integral trapezioid = " << integral2;
		integral1 = integralT(a, b, eps, n * 2);
		R = (integral1 - integral2) * (pow(0.5, 2) - 1);
		cout << endl << "error = " << R;
	}
	if (choice == 2) {
		integral2 = integralS(a, b, eps, n);
		cout << "Integral Simpson = " << integral2;
		integral1 = integralS(a, b, eps, n * 2);
		R = (integral1 - integral2) * (pow(0.5, 4) - 1);
		cout << endl << "error = " << R;
	}
	if (choice != 2 && choice != 1) {
		cout << "wrong input";
	}
	double h, k, a1, b1, c, d;
	int nx = 21, ny = 21; // количество шагов
	a1 = 0, b1 = 2, c = 0, d = 1;
	cout << "Double integral =" << doubleIntegral(nx, ny, a1, b1, c, d);
}

double function(double x) {
	return ((1 + pow(x, 2)) / (1 + pow(x, 3)));
}
double function2(double x, double y) {
	return (pow(x, 2) + y);
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
	} while (abs(Integr1 - Integr2) >= 3 * eps);
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
	return Integr1;
}
double doubleIntegral(double nx, double ny, double lx, double ux, double ly, double uy) {
	double h, k;
    // в z храним таблицу значений
    // ax[] хранит интегральные значения y для всех рассмотренных х
    double z[50][50], ax[50], answer;
    h = (ux - lx) / (nx-1);// рассчет шага
	k = (uy - ly) / (ny-1);
    for (int i = 0; i < nx; ++i) { // рассчет значений функции для таблицы
        for (int j = 0; j < ny; ++j) {
            z[i][j] = function2(lx + i * h, ly + j * k);
        }
    } 
    // вычисление интеграла относительно y в каждой точке для x *
    for (int i = 0; i < nx; ++i) {
        ax[i] = 0;
        for (int j = 0; j < ny; ++j) {
            if (j == 0 || j == ny - 1) 
                ax[i] += z[i][j];
            else if (j % 2 == 0)
                ax[i] += 2 * z[i][j]; //четные * 2
            else
                ax[i] += 4 * z[i][j]; // нечетные
        }
        ax[i] *= (k / 3);
    }
    answer = 0;
    // Вычисление итогового значения интеграла используя значения полученные раньше
    for (int i = 0; i < nx; ++i) {
        if (i == 0 || i == nx - 1)
            answer += ax[i]; //первая и последняя
        else if (i % 2 == 0)  //если четные номера, то * 2, если нет, то *4
            answer += 2 * ax[i];
        else
            answer += 4 * ax[i];
    }
    answer *= (h / 3);
    return answer;
}
