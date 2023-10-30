//решить систему уравнени1 методом ньютона
#include <iostream>
#include <math.h>
#include<iomanip>
using namespace std;
double function1(double x1, double x2);
double function2(double x1, double x2);
double f1derivative11(double x1, double x2, double dif);
double f1derivative12(double x1, double x2, double dif);
double f2derivative21(double x1, double x2, double dif);
double f2derivative22(double x1, double x2, double dif);
void residualVector(double* vNev, double x1, double x2);
double* gauss(double** matrix, double* y, const int n, float accur);
void outputMatrix(double** matrix, double* y, const int n);

void main() {
	double x1 = 1, x2 = 0;
    int k=1, N=100;
    double accur = 1e-9, dif = 0.001; 
	const int n = 2;
	double* vNev, *xVect, *delta, *solution, *vNev2;
	vNev = new double[n];
    vNev2 = new double[n];
    delta = new double[n];
    xVect = new double[n];
    solution = new double[n];
    solution[0] = x1;
    solution[1] = x2;
	double** a = new double* [n];
	for (int i = 0; i < n; i++) { //создаем копии
		a[i] = new double[n];
	}
    double** a2 = new double* [n];
    for (int i = 0; i < n; i++) {
        a2[i] = new double[n];
    }
    do{
        residualVector(vNev, x1, x2); //заполнение
	    a[0][0] = f1derivative11(x1, x2, dif);
	    a[1][0] = f1derivative12(x1, x2, dif);
	    a[0][1] = f2derivative21(x1, x2, dif);
	    a[1][1] = f2derivative22(x1, x2, dif);
        for (int i = 0; i < n; i++) { //копируем матрицу
            for (int j = 0; j < n; j++) {
                a2[i][j] = a[i][j];
            }
        }
        for (int i = 0; i < n; i++) {//копируем вектор невязки
            vNev2[i] = vNev[i];
        }
        outputMatrix(a2, vNev2, n);
        delta = gauss(a2, vNev, n, accur); //дельта - решение уравнения
        for (int i = 0; i < n; i++)
            solution[i] += delta[i]; //к решению(первое приближение) прибавляется дельта
        for (int i = 0; i < n; i++) {
            cout << "delta" << i+1 << " = "<< delta[i] << " " << endl;
        }
        double max1 = 0;
        double max2 = 0;
        for (int i = 0; i < n; i++)  // максимально маленькие дельта
        {
            if (abs(vNev[i]) > max1) {
                max1 = abs(vNev[i]);
            }
            if (abs(solution[i]) < 1)
            {
                if (abs(delta[i]) > max2)
                    max2 = abs(delta[i]);
            }
            if (abs(delta[i] >= 1))
            {
                if (abs(delta[i] / solution[i]) > max2)
                    max2 = abs(delta[i]);
            }
        }

            delta[0] = max1;
            delta[1] = max2;
            cout << endl;
            x1 = solution[0];
            x2 = solution[1];
            k++;
            if (k >= N) 
            {
            break;
            }

    } while (delta[0] > accur || delta[1] > accur); //точность

    x1 = solution[0];
    x2 = solution[1];
    cout << "____________________________________________________________" << endl;
    for (int i = 0; i < n; i++)
        cout << solution[i] << endl; //удаление
    delete[] vNev;  
    delete[] vNev2;
    delete[] xVect;
    delete[] delta;
    delete[] solution;
    for (int i = 0; i < n; i++) {
        delete[] a[i];
    }
    for (int i = 0; i < n; i++) {
        delete[] a2[i];
    }
}

double function1(double x1, double x2) {
	return (sin(x1) - x2 - 1.32);
}
double function2(double x1, double x2) {
	return (cos(x2) - x1 + 0.85);
}
double f1derivative11(double x1, double x2, double dif) {
    cout << ((function1(x1 + dif, x2) - function1(x1, x2)) / (dif));
	return((function1(x1 + dif, x2) - function1(x1, x2)) / (dif));
}
double f1derivative12(double x1, double x2, double dif) {
    cout << ((function1(x1, x2 + dif) - function1(x1, x2)) / (dif));
	return((function1(x1, x2 + dif) - function1(x1, x2)) / (dif));
}
double f2derivative21(double x1, double x2, double dif) {
    cout << ((function2(x1 + dif, x2) - function2(x1, x2)) / (dif));
	return((function2(x1 + dif, x2) - function2(x1, x2)) / (dif));
}
double f2derivative22(double x1, double x2, double dif) {
    cout << ((function2(x1, x2 + dif) - function2(x1, x2)) / (dif));
	return((function2(x1, x2 + dif) - function2(x1, x2)) / (dif));
}
//double f1derivative11(double x1, double x2) {
//	return(cos(x1));
//}
//double f1derivative12(double x1, double x2) {
//	return(-1);
//}
//double f2derivative21(double x1, double x2) {
//	return(-1);
//}
//double f2derivative22(double x1, double x2) {
//	return(-sin(x2));
//}
void residualVector(double* vNev, double x1, double x2)
{
	vNev[0] = -function1(x1, x2);
	vNev[1] = -function2(x1, x2);
}
void outputMatrix(double** matrix, double* y, const int n) {
    cout << endl << "matrix: " << endl;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {

            cout << matrix[j][i] << setw(10);

        }
        cout << "|" << y[j] << " " << endl;
    }
}
double* gauss(double** matrix, double* y, const int n, float accur)
{
    double* xVect, max, roundH;
    int k, index;
    const double eps = 0.00001;  // точность
    xVect = new double[n];
    k = 0;
    while (k < n)
    {
        // Поиск максимума
        max = abs(matrix[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(matrix[i][k]) > max)
            {
                max = abs(matrix[i][k]);
                index = i;
            }
        }
        for (int j = 0; j < n; j++)
        {
            swap(matrix[k][j], matrix[index][j]);
        }
        swap(y[k], y[index]);
        if (k > 0)
            outputMatrix(matrix, y, n);
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = matrix[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++) {
                matrix[i][j] = matrix[i][j] / temp;
            }
            y[i] = y[i] / temp;           
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            y[i] = y[i] - y[k];
        }
        if (k == n - 1) {
            outputMatrix(matrix, y, n);
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        xVect[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - matrix[i][k] * xVect[k];
    }
    return xVect;
}
