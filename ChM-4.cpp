#include <iostream>
#include <math.h>
#include<iomanip>
using namespace std;
double sumT1(double* t, int amount);
double sumT2(double* t, int amount);
double sumU1(double* u, int amount);
double sumU2(double* t, double* u, int amount);
double* gauss(double** matrix, double* y, const int n);
void outputMatrix(double** matrix, double* y, const int n);
void makeMatrix(double** matrix, double* vNev, double* t, double* u, int amount);

void main() {
    //	n=c*F^(-1/e) e, c - ? n = c/êîðåíü ñòåïåíè å èç F //îòìåíà // âàðèàíò 5:    u=a*exp(-Bt) a è B - ? --> ln(u)=ln(a)-Bt
    double* tString, * uString, * vNev, * xVect, *answerValue;
    int n = 7, nMatr = 2;
    answerValue = new double[nMatr];
    vNev = new double[nMatr];
    xVect = new double[nMatr];
    tString = new double[n];
    tString[0] = 60; tString[1] = 70; tString[2] = 80; tString[3] = 90; tString[4] = 100; tString[5] = 110; tString[6] = 120;
    uString = new double[n];
    uString[0] = 0.0148; uString[1] = 0.0124; uString[2] = 0.0102; uString[3] = 0.0085; uString[4] = 0.0071; uString[5] = 0.0059; uString[6] = 0.0051;
    double** matrix = new double* [nMatr];
    for (int i = 0; i < nMatr; i++) { //ñîçäàåì êîïèè
        matrix[i] = new double[nMatr];
    }

    makeMatrix(matrix, vNev, tString, uString, n);
    xVect = gauss(matrix, vNev, nMatr);
    cout << xVect[0] << " " << xVect[1]<< endl; //âîçâðàùàåò ln(a) è -B
    double sumS=0;
    answerValue[0] = exp(xVect[0]);
    answerValue[1] = (-1)*xVect[1];
    cout << "a = " << answerValue[0] << "; b = " << answerValue[1]<<endl;
    
    for (int i = 0; i < n; ++i) {
        sumS += pow((log(uString[i]) - xVect[0] - xVect[1] * tString[i]), 2);
    }
    double S2 = (1 / (n - 2)) * sumS;
    cout << pow(S2, 1 / 2);
}
double sumT1(double* t, int amount) {
    double sum = 0;
    for (int i = 0; i < amount; ++i) {
        sum += t[i];
    }
    return sum;
}
double sumT2(double* t, int amount) {
    double sum = 0;
    for (int i = 0; i < amount; ++i) {
        sum += pow(t[i], 2);
    }
    return sum;
}
double sumU1(double* u, int amount) {
    double sum = 0;
    for (int i = 0; i < amount; ++i) {
        sum += (log(u[i]));
    }
    return sum;
}
double sumU2(double* t, double* u, int amount) {
    double sum = 0;
    for (int i = 0; i < amount; ++i) {
        sum += t[i] * (log(u[i]));
    }
    return sum;
}
double* gauss(double** matrix, double* y, const int n)
{
    double* xVect, max, roundH;
    int k, index;
    const double eps = 0.00001;  // òî÷íîñòü
    xVect = new double[n];
    k = 0;
    while (k < n)
    {
        // Ïîèñê ìàêñèìóìà
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
        // Íîðìàëèçàöèÿ óðàâíåíèé
        for (int i = k; i < n; i++)
        {
            double temp = matrix[i][k];
            if (abs(temp) < eps) continue; // äëÿ íóëåâîãî êîýôôèöèåíòà ïðîïóñòèòü
            for (int j = 0; j < n; j++) {
                matrix[i][j] = matrix[i][j] / temp;
            }
            y[i] = y[i] / temp;
            if (i == k)  continue; // óðàâíåíèå íå âû÷èòàòü ñàìî èç ñåáÿ
            for (int j = 0; j < n; j++)
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            y[i] = y[i] - y[k];
        }
        if (k == n - 1) {
            outputMatrix(matrix, y, n);
        }
        k++;
    }
    // îáðàòíàÿ ïîäñòàíîâêà
    for (k = n - 1; k >= 0; k--)
    {
        xVect[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - matrix[i][k] * xVect[k];
    }
    return xVect;
}
void makeMatrix(double** matrix, double* vNev, double* t, double* u, int amount) {
    matrix[0][0] = 7; matrix[1][0] = sumT1(t, amount);
    matrix[0][1] = sumT1(t, amount); matrix[1][1] = sumT2(t, amount);
    vNev[0] = sumU1(u, amount);  vNev[1] = sumU2(t, u, amount);
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
