#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace std;

// вывести матрицу
void printingA(int n, int t, double** z)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < t; j++)
		{
			std::cout << setw(3) << z[i][j];
		}
		std::cout << endl;
	}
	std::cout << endl;
}
//вывести вектор
void printingV(int n, double* z)
{
	int i;
	for (i = 0; i < n; i++)
	{
		
		std::cout <</* setw(7) <<*/ z[i] << endl;
		
	}
	std::cout << endl;
}
// транспонировать матрицу
double** transp(double** a, int n, int m)
{

	int i, j;
	double** b;
	b = new double* [m];
	for (i = 0; i < m; i++)
	{
		b[i] = new double[n];
	}

	//написать транспонирование(матрица а расширенная!!!!!!!!!!!!!!!)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			b[i][j] = a[j][i];
	}
	//printingA(m, n, b);
	return b;
}

/*генерация хорошо обусловленных матриц с правой частью
для единичного решения*/
/*матрица положительна определена и симметрична*/
double** generateMatrix(int n, double*& f)
{
	srand(time(NULL));
	double** matrix = new double* [n],
		** L = new double* [n];
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		L[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) L[i][j] = 0.0;
			if (i == j) L[i][j] =  1.0137*(rand()%50);
			if (i > j) L[i][j] = 1.0137 * (rand()%100)-50;
		}
	}
	
	//printingA(n, n, L);
	double** LT = transp(L, n, n);
	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
		{
			matrix[i][j] = 0.0;
			for (int k = 0; k < n; k++)
			{
				matrix[i][j] += LT[j][k] * L[k][i];
			}
		}
	}
	//printingA(n, n, matrix);
	for (int i = 0; i < n; i++)
	{
		f[i] = 0.0;
		for (int j = 0; j < n; j++)
		{
			f[i] += matrix[i][j];
		}
	}
	//printingV(n, f);
	return matrix;
}

/*генерация плохообусловленной матрицы размера n и порядка k
с правой частью для единичного решения*/
double** generateIllCondition(int n, int k, double*& f)
{
	srand(time(NULL));
	/*создаем нижнетреугольную
	и верхнетреугольную матрицу*/
	double** matrix = new double* [n],
		** L = new double* [n];
	for (int i = 0; i < n; i++)
	{
		matrix[i] = new double[n];
		L[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) L[i][j] = 0.0;
			if (i == j) L[i][j] = ((double)rand() / RAND_MAX) * 20000;
			if (i > j) L[i][j] = ((double)rand() / RAND_MAX) * 20000 - 10000;
		}
	}
	double** LT = transp(L, n, n);
	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
		{
			matrix[i][j] = 0.0;
			for (int k = 0; k < n; k++)
			{
				matrix[i][j] += LT[j][k] * L[k][i];
			}
		}


	}
	for (int i = 0; i < n; i++)
	{
		matrix[i][i] *= pow(10, -k);
	}
	for (int i = 0; i < n; i++)
	{
		f[i] = 0.0;
		for (int j = 0; j < n; j++)
		{
			f[i] += matrix[i][j];
		}
	}
	return matrix;

}

double* metodSquareRoot(double** a, int N, int S, double* f, double* b, double*& r, double& p, bool& flag)
{

	double** raba = new double* [S];
	double* rabf = new double[S];
	for (int i = 0; i < S; i++)
	{
		raba[i] = new double[S];
	}
	double** tempa = new double* [N];
	double* tempf = new double[N];
	for (int i = 0; i < N; i++)
	{
		tempa[i] = new double[S];
		tempf[i] = f[i];
	}
	if (N != S)
	{
		//матрицу а портить нельзя
		//и исходную правую часть тоже

		
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < S; j++)
			{
				tempa[i][j] = a[i][j];

			}
		}
		double** aT = transp(a, N, S);
		for (int i = 0; i < S; i++)
		{
			for (int j = 0; j < N; j++)
			{
				tempa[j][i] = tempa[j][i] * b[j];
				if (i == 0)tempf[j] *= b[j];
			}
		}
		//printingA(N, S, tempa);
		for (int i = 0; i < S; i++)
		{
			for (int j = 0; j < S; j++)
			{
				raba[i][j] = 0.0;
				for (int k = 0; k < N; k++)
				{
					raba[i][j] += aT[j][k] * tempa[k][i];
				}
			}
		}
		for (int i = 0; i < S; i++)
		{
			double sum = 0;
			for (int j = 0; j < N; j++)
			{
				sum += aT[i][j] * tempf[j];
			}
			rabf[i] = sum;
		}

		//printingA(S, S, raba);
		//printingV(S, rabf);

	}
	else
	{
		for (int i = 0; i < S; i++)
		{
			for (int j = 0; j < S; j++)
			{
				raba[i][j] = a[i][j];
				if (i == 0) rabf[j] = f[j];
			}
		}
	}

	double** T = new double* [S];
	for (int i = 0; i < S; i++)
	{
		T[i] = new double[S];
	}
	for (int i = 0; i < S; i++)
	{
		for (int j = 0; j < S; j++)
		{
			T[i][j] = 0;
		}
	}
	for (int i = 0; i < S; i++)
	{
		double sum = 0.0;
		for (int k = 0; k < i; k++)
		{
			sum += T[k][i] * T[k][i];
		}
		T[i][i] = sqrt(raba[i][i] - sum);
		for (int j = i ; j < S; j++)
		{
			double sum = 0.0;
			for (int k = 0; k < i; k++)
			{
				sum += T[k][i] * T[k][j];
			}
			T[i][j] = (raba[i][j] - sum) / T[i][i];
		}
	}
	//printingA(S, S, T);
	
	//double** TT = transp(T, S, S);
	//printingA(S, S, TT);
	 //std::cout << "tt " << T[S-1][S-1] <<endl;
	if (!_isnan(T[S-1][S-1]) && T[S-1][S-1]!=0 /*&& abs(T[S-1][S-1])>=1e-7*/) {
		
		double* y = new double[S];
		double* x = new double[S];
		for (int i = 0; i < S; i++)
		{
			for (int j = 0; j < i; j++)
			{
				rabf[i] -= T[j][i] * y[j];
			}
			y[i] = rabf[i] / T[i][i];
			
		}
		//std::cout << y[S - 1] <<" "<<rabf[S-1]<<" "<<T[S-1][S-1]<< endl;
		//printingV(S, y);
		for (int i = S - 1; i >= 0; i--)
		{
			for (int j = i + 1; j < S; j++)
			{
				y[i] -= T[i][j] * x[j];
			}
			x[i] = y[i] / T[i][i];
		}
		//printingV(S, x);
		for (int i = 0; i < N; i++)
		{
			r[i] = 0;
			for (int j = 0; j < S; j++)
			{
				r[i] += a[i][j] * x[j];
			}
			r[i] = abs(r[i] - f[i]);
		}

		p = 0;
		for (int i = 0; i < N; i++)
		{
			p += r[i] * r[i];
		}
		p = sqrt(p);
		return x;
	}
	else 
	{
		flag = false;
		return NULL;
	}

	
}

double test(int numberTest, int n, int k)
{
	double maxInac = 0.0;
	double inaccuracy = 0.0;
	for (int i = 0; i < numberTest;)
	{
		double* f = new double[n];
		double* b = new double[n];
		double* r = new double[n];
		double* x = new double[n];
		double p = 0.0;
		bool flag = true;
		if (k != 0)
		{
			double** a = generateIllCondition(n, k, f);
			 x = metodSquareRoot(a, n, n, f, b, r, p,flag);
			 
			delete[]a;
		}
		else
		{
			double** a = generateMatrix(n,f);
			 x = metodSquareRoot(a, n, n, f, b, r, p,flag);
			 

			
		}
		delete []  b;
		delete [] r;
		delete [] f;
		//printingV(n, x);
		//std::cout << sizeof(x) / sizeof(double) << endl;
		if (flag)
		{
			maxInac = abs(x[0] - 1);
			for (int i = 1; i < n; i++)
			{
				maxInac = abs(x[i] - 1) > maxInac ? abs(x[i] - 1) : maxInac;
			}
			//if (maxInac > 1) printingV(n, x);
			delete[] x;
			inaccuracy += maxInac;
			i++;
		}
	}
	
	
	return inaccuracy / numberTest;
}
int main()
{
	setlocale(LC_ALL, "RUS");
	std::cout << "Численный методы.Лабораторная работа №2" << endl << "Выполнил Мазуров А.А." << endl << "Проверила Шабунина З.А. \n";
	std::cout << "Тестовые примеры для переопределенных систем уравнений" << endl;
	int n = 3, s = 2;
	bool flag = true;
	double** a = new double* [n];
	for (int i = 0; i < n; i++)
	{
		a[i] = new double[s];
	}
	a[0][0] = 1.0;
	a[0][1] = 0.0;
	a[1][0] = 0.0;
	a[1][1] = 1.0;
	a[2][0] = 2.0;
	a[2][1] = 3.0;
	double f[3] = { 0,0,8 };
	std::cout <<"Тестовая матрица " <<endl;
	printingA(n, s, a);
	std::cout<<"Вектор правой части"<<endl;
	printingV(n, f);
	double* r1 = new double[n],
		* r2 = new double[n],
		* r3 = new double[n];
	double p1 = 0.0, p2 = 0.0, p3 = 0.0;
	
	double b1[3] = { 1,1,1 }, b2[3] = { 2,2,1 },
		b3[3] = { 1,1,2 };
	double* x1 = metodSquareRoot(a, n, s, f, b1, r1, p1,flag);
	std::cout << "Весовые коэффициенты" << endl;
	printingV(n, b1);
	std::cout << "Вектор решения " << endl;
	printingV(s,x1);
	std::cout << "Вектор невязки " << endl;
	printingV(n, r1);
	std::cout << "Норма вектора невязки " << p1 << endl;
	std::cout << endl;
	double* x2 = metodSquareRoot(a, n, s, f, b2, r2, p2,flag);
	std::cout << "Весовые коэффициенты" << endl;
	printingV(n, b2);
	std::cout << "Вектор решения " << endl;
	printingV(s, x2);
	std::cout << "Вектор невязки " << endl;
	printingV(n, r2);
	std::cout << "Норма вектора невязки " <<p2<< endl;
	std::cout << endl;
	double* x3 = metodSquareRoot(a, n, s, f, b3, r3, p3,flag);
	std::cout << "Весовые коэффициенты" << endl;
	printingV(n, b3);
	std::cout << "Вектор решения " << endl;
	printingV(s, x3);
	std::cout << "Вектор невязки " << endl;
	printingV(n, r3);
	std::cout << "Норма вектора невязки " << p3 << endl;
	std::cout << endl;

	std::cout << "Тестирование для хорошо обусловленных матриц" << endl;
	srand(time(NULL));
	std::cout << (double)(rand() % 500) << " " << (double)(rand() % 1000) - 500 << endl;
	for (int i = 20; i < 150; i += 20)
	{
		std::cout<<"Размерность "<<i<<" "<<endl;
		double inaccuracy=test(3, i, 0);
		std::cout << "Средняя относительная погрешность " << inaccuracy << endl;
	}
	/*std::cout << endl;
	std::cout << "Тестирование для плохо обусловленных матриц" << endl;
	for (int i = 10; i < 30; i += 10)
	{
		std::cout << "Размерность матрицы " <<i<< endl;
		for (int j = 2; j < 8; j += 2)
		{
			std::cout << "Порядок k " <<j<<" ";
			double inaccuracy = test(5, i, j);
			std::cout << "Средняя относительная погрешность " << inaccuracy << endl;

		}
	}*/



















	//printing(1, s, x3);
	//printing(1, n, f);
	/*if (!flag)
	{
		cout << "Матрица не симметричная, необходимо домножить на транспонированную к ней" << endl;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				at[j][i] = a[i][j];
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				rab[i][j] = 0;
				b1[i] = 0;
				for (k = 0; k < n; k++)
				{
					rab[i][j] += at[i][k] * a[k][j];
					b1[i] += at[i][k] * b[k];
				}
			}
		}
	}
	else
	{
		cout << "Матрица симметричная" << endl;
		for (i = 0; i < n; i++)
		{
			b1[i] = b[i];
			for (j = 0; j < n; j++)
				rab[i][j] = a[i][j];
		}
	}
	cout << "Рабочая матрица:" << endl;
	printing(n, n, rab);
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < (i + 1); k++)
		{
			double sum = 0;
			for (int j = 0; j < k; j++)
				sum += s[i][j] * s[k][j] * d[j][j];
			if (i == k)
			{
				s[i][k] = sqrt(rab[i][i] - sum);
				d[i][k] = sign(rab[i][i] - sum);
			}
			else
				s[i][k] = (1.0 / s[k][k] * (rab[i][k] - sum));
		}
	}
	s = transp(s, n, n);
	cout << "Mатрица S:" << endl;
	printing(n, n, s);
	cout << "Вектор Y: (";
	for (i = 0; i < n; i++)
	{

		double sum = 0;
		for (int k = 0; k <= i - 1; k++)
			sum += y[k] * s[k][i];
		y[i] = (b1[i] - sum) / s[i][i];
		cout << y[i];
		if (i != n - 1)
			cout << ", ";
		else
			cout << ")";
	}
	cout << endl;
	cout << "Вектор X: (";
	for (i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int k = i + 1; k <= n - 1; k++)
			sum += s[i][k] * x[k];
		x[i] = (y[i] - sum) / s[i][i];
	}
	for (i = 0; i < n; i++)
	{
		opred *= s[i][i] * s[i][i];
		cout << x[i];
		if (i != n - 1)
			cout << ", ";
		else
			cout << ")";
	}
	cout << endl;
	cout << "Определитель исходной матрицы: " << opred << endl;
	cout << "Вектор невязки: (";
	for (i = 0; i < n; i++)
	{
		nev[i] = -b[i];
		for (k = 0; k < n; k++)
			nev[i] += a[i][k] * x[k];
		cout << nev[i];
		if (i != n - 1)
			cout << ", ";
		else
			cout << ")";
	}
	cout << endl;*/
	system("pause");
}