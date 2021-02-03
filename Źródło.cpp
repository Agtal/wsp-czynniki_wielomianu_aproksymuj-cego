//Agnieszka Talik
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <mkl_lapack.h>

using namespace std;

#define size 8 // wielkosc wczytanej tablicy danych
#define bsize 2 // stopien wielomianu

double x[size], y[size];
int w = 1; // waga
int n = 0;

void make_tab(double a, double b)
{
	x[n] = a;
	y[n] = b;

	n++;
}

void matrix(double x[], double y[], double gkj[bsize][bsize], double Fk[])
{
	double aj=0, ak=0;
	for (int k = 0; k < bsize; k++)
		{
		for (int j = 0; j < bsize; j++)
			{
				for (int i = 0; i < size; i++)
				{
					ak = pow(x[i], k);
					aj = pow(x[i], j);

					gkj[k][j] += ak * aj * w;
				}
				cout << "g" << k << j << " = " << gkj[k][j] << setw(10);
			}
				cout << setw(0) << endl;
		}
}

void arrayF(double x[], double y[], double Fk[])
{
	cout << endl << endl ;
	double ak = 0;
	for (int k = 0; k < bsize; k++)
	{
		for (int i = 0; i < size; i++)
		{
			ak = pow(x[i], k);
			Fk[k] += ak * y[i] * w;
		}
		cout << "F[" << k << "] = " << Fk[k] << endl;
	}
	cout << endl << endl;
}

void G(double* Fk, double* x)
{
	cout << endl << endl << "Funkcja aprokcymujaca: G(x) = " << Fk[0] << " * x^0 + " << Fk[1] << " * x^1" << endl << endl;
	double g[size] = {};
	for (int i = 0; i < size; i++)
	{
		g[i] = (Fk[0] * pow(x[i], 0)) + (Fk[1] * pow(x[i], 1));
		cout << "G(x" << i << ") = " << g[i] << endl;
	}
}

int main() {

	ifstream plik;
	plik.open("tabela.txt");

	while (true)
	{
		double a, b;
		plik >> a >> b;
		if (!plik.fail())
			make_tab(a, b);
		else
			break;
	}
	plik.close();
	cout << "n" << setw(5) << "x" << setw(5) << "y" << endl;


	for (int i = 0; i < size; i++)
	{
		cout << i << setw(5) << x[i] << setw(5) << y[i] << endl;
	}

	cout << endl << "Liczba wezlow = "<< n << endl <<
		"Stopien wielomianu = " << bsize - 1 << endl << endl;

	double gkj[bsize][bsize] = {};
	double Fk[bsize] = {};
	matrix(x, y, gkj, Fk);
	arrayF(x, y, Fk);
		
	int info;
	// lda = ldb = N = bsize -  Liczba równañ liniowych, czyli rz¹d macierzy A. N> = 0.
	int N = bsize;
	int NRHS = 1; // liczba kolumn macierzy

	dposv("U", &N, &NRHS, *gkj, &N, Fk, &N, &info);

	cout << "Wspolczynniki wielomianu: " << endl;
	for (int k = 0; k < bsize; k++)
	{
		cout << "a[" << k << "] = " << Fk[k] << endl;
	}

	G(Fk, x);
}