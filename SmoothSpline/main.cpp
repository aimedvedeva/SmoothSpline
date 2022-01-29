#include <iostream>
#include <fstream>

using namespace std;
#define pi 3.1415926536

ofstream fileOut;

double ContinuousFun(double value) {
	return sinl(value*pi);
}

double PiecewiseFun(double value) {
	if (value < -1.25) {
		return 3 * sqrt(abs(value - 1)) - 3.25;
	}
	else if (-1.25 <= value && value < 1.25) {
		return -value;
	}
	else if (value >= 1.25) {
		return 3.35 - 3 * sqrt(value + 1);
	}
}

//Uniform grid
double *GetEvenGrid(int n, double a, double b, double *delta) {
	double *grid = new double[n + 1];
	*delta = (b - a) / (double)n;

	grid[0] = a;
	for (int i = 1; i < n + 1; i++) {
		grid[i] = grid[i - 1] + *delta;
	}

	return grid;
}

//Grid function
double *GetFunValuesGrid(double *grid, double (*Fun) (double), int n) {
	double *funGrid = new double[n];

	for (int i = 0; i < n; i++) {
		funGrid[i] = Fun(grid[i]);
	}
	return funGrid;
}

//Functions for working with matrices
void MultMatrixByNum(double *matrix, int numOfColumns, int numOfRows, double num) {
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrix[i * numOfColumns + j] *= num;
		}
	}
}

void InitMatixByZero(double *matrix, int width, int height) {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			matrix[i * width + j] = 0;
		}
	}
}

double *MultMatrixAB(double *matrixA, int numOfVariablesA, int numOfEquationsA, double *matrixB, int numOfVariablesB, int numOfEquationsB) {
	if (numOfVariablesA != numOfEquationsB) {
		return NULL;
	}

	double *newMatrix = new double[numOfEquationsA * numOfVariablesB];

	for (int i = 0; i < numOfEquationsA; i++) {
		for (int j = 0; j < numOfVariablesB; j++) {
			double sum = 0;
			for (int k = 0; k < numOfVariablesA; k++) {
				sum += matrixA[i * numOfVariablesA + k] * matrixB[k * numOfVariablesB + j];
			}
			newMatrix[i * numOfVariablesB + j] = sum;
		}
	}

	return newMatrix;
}

double *TransposeMatrix(double *matrix, int numOfColumns, int numOfRows)
{
	double *newMatrix = new double[numOfColumns * numOfRows];
	int counter = 0;
	for (int j = 0; j < numOfColumns; j++) {
		for (int i = 0; i < numOfRows; i++) {
			newMatrix[counter] = matrix[i * numOfColumns + j];
			counter++;
		}
	}

	return newMatrix;
}

void InverseDiagMatrix(double *matrix, int width, int height) {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (i == j) {
				matrix[i * width + j] = 1 / matrix[i * width + j];
			}
		}
	}
}

double *AddMatrix(double *matrixA, double *matrixB, int numOfColumns, int numOfRows) {
	double *sumMatrix = new double[numOfColumns * numOfRows];
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			sumMatrix[i * numOfColumns + j] = matrixA[i * numOfColumns + j] + matrixB[i * numOfColumns + j];
		}
	}
	return sumMatrix;
}

//The functions are named in terms of the theoretical SLAE: (A + H * P^(-1) * H^(T)) * m = H * y;   A + H * P^(-1) * H^(T) = R;   H * y = B;  => R * m = B -> sweep method
double *GetAMatix(int n, double delta) {
	int width = n - 1;
	int height = n - 1;
	double *A = new double[width * height];
	InitMatixByZero(A, width, height);
	int rowCounter = 0;

	for (int i = 0; i < height; i++) {
	
		if (i != 0) {
			A[i * width + rowCounter - 1] = delta / 6;
		}

		A[i * width + rowCounter] = 2 * delta / 3;

		if (i != width - 1) {
			A[i * width + rowCounter + 1] = delta / 6;
		}
	    rowCounter++;
	}
	return A;
}

double *GetHMatrix(int n, double delta) {
	int width = n + 1;
	int height = n - 1;
	double *H = new double[width * height];
	InitMatixByZero(H, width, height);
	int rowCounter = 0;

	for (int i = 0; i < height; i++) {
		
		H[i * width + rowCounter] = 1 / delta;
		H[i * width + rowCounter + 1] = - 2 / delta;
		if (i * width + rowCounter + 2 < width * height) {
			H[i * width + rowCounter + 2] = 1 / delta;
		}
     	
		rowCounter++;
	}
	return H;
}

double *GetPMatrix(int n, double value) {
	int width = n + 1;
	int height = n + 1;
	double *P = new double[width * height];
	InitMatixByZero(P, width, height);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (i == j) {
				P[i * width + j] = value;
			}
		}
	}
	return P;
}

double *GetRMatrix(double *H, int Hwidth, int Hheight, double *A, int Awidth, int Aheight, double alpha) {
	double *transposedH = TransposeMatrix(H, Hwidth, Hheight);
	double *HtransposedH = MultMatrixAB(H, Hwidth, Hheight, transposedH, Hheight, Hwidth);
	delete[] transposedH;

	MultMatrixByNum(HtransposedH, Hheight, Hheight, alpha / (1 - alpha));
	
	double *sumAHH = AddMatrix(A, HtransposedH, Awidth, Aheight);
	delete[] HtransposedH;

	return sumAHH;
}

double *GetBMatrix(double *H, int Hwidth, int Hheight, double *funValuesVector, int lenght) {
	double *HFunValVec = MultMatrixAB(H, Hwidth, Hheight, funValuesVector, 1, lenght);
	return HFunValVec;
}

void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	fileOut.precision(10);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			fileOut << matrix[i * numOfVariables + j] << "     ";
		}
		fileOut << endl;
	}
	fileOut << endl;
}

void PrintNum(double num){
	fileOut << num << endl;
}


int main(void) {
	int numOfDots = 7;
	double alpha = 0.1;
	double delta;
	double *grid = GetEvenGrid(numOfDots, 0, 2, &delta);
	double *funValGrid = GetFunValuesGrid(grid, ContinuousFun, numOfDots + 1);
	double *H = GetHMatrix(numOfDots, delta);
	
	ofstream fileH;
	fileH.open("H.txt");
	fileH.precision(10);

	for (int i = 0; i < numOfDots - 1; i++) {
		for (int j = 0; j < numOfDots + 1; j++) {
			fileH << H[i * (numOfDots + 1) + j] << "     ";
		}
		fileH << endl;
	}
	fileH << endl;
	fileH.close();
	
	double *A = GetAMatix(numOfDots, delta);
	double *R = GetRMatrix(H, numOfDots + 1, numOfDots - 1, A, numOfDots - 1, numOfDots - 1, alpha);                                                                         
	double *B = GetBMatrix(H, numOfDots + 1, numOfDots - 1, funValGrid, numOfDots + 1);

	fileOut.open("output.txt");
	PrintNum(numOfDots);
	PrintNum(delta);
	PrintNum(alpha);
	PrintMatrix(funValGrid, 1, numOfDots + 1);
	PrintMatrix(grid, 1, numOfDots + 1);
	PrintMatrix(R, numOfDots - 1, numOfDots - 1);
	PrintMatrix(B, 1, numOfDots - 1);
	//PrintMatrix(A, numOfDots - 1, numOfDots - 1);
	fileOut.close();

	delete[] B;
	delete[] R;
	delete[] A;
	delete[] H;
	delete[] funValGrid;
	delete[] grid;
}

