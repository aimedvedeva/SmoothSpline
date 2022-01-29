#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <crtdbg.h>

using namespace std;


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

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
    }
	sum = sqrt(sum);
	return sum;
}

double *TransposeMatrix(double *matrix, int numOfColumns, int numOfRows)
{
	double *newMatrix = new double[numOfColumns * numOfRows];

	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			newMatrix[j * numOfColumns + i] = matrix[i * numOfColumns + j];
		}
	}

	return newMatrix;
}
/*
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
}*/

void SubtractMatrix(double *matrixA, double *matrixB, int numOfColumns, int numOfRows) {
	
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrixA[i * numOfColumns + j] -= matrixB[i * numOfColumns + j];
		}
	}
}

void CopyMatrix(double *to, double *from, int numOfRows, int numOfColumns, int fromColumn, int fromRow) {
	for (int i = fromRow; i < numOfRows; i++) {
		for (int j = fromColumn; j < numOfColumns; j++) {
			to[i * numOfColumns + j] = from[i * numOfColumns + j];
		}
	}
}

double *Normal(double *matrixA, int numOfColumns, int numOfEquations, int step) {
	double *normal = new double[numOfEquations - step];

	//копируем нужный нам столбец
	for (int i = step; i < numOfEquations; i++) {
		normal[i - step] = matrixA[i * numOfColumns + step];
	}

	double factor = EuclideanNorm(normal, numOfEquations - step);

	if (normal[0] >= 0) {
		factor = -factor;
	}

	normal[0] -= factor;//имитации вычитания первого вектора евстественного базина * factor

	return normal;
}

void InitByZero(double *data, int size) {
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}

double *GetIdentityMatrix(int numOfColumns, int numOfRows) {
	double *identityMatrix = new double[numOfColumns * numOfRows];
    
	
	InitByZero(identityMatrix, numOfColumns * numOfRows);
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			if (i == j) {
				identityMatrix[i * numOfColumns + j] = 1;
			}
		}
	}
	return identityMatrix;
}

void MultMatrixByNum(double *matrix, int numOfColumns, int numOfRows, double num) {
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrix[i * numOfColumns + j] *= num;
		}
	}
}

double *ReflectionMatrix(double *normal, int lenght) {
	double norm = EuclideanNorm(normal, lenght);

	double *reflectionMatrix = GetIdentityMatrix(lenght, lenght);

	double *mult = MultMatrixAB(normal, 1, lenght, normal, lenght, 1);//w * wT

	MultMatrixByNum(mult, lenght, lenght, (2 / powl(norm, 2)));//2/norm^2 * (w * wT)

	SubtractMatrix(reflectionMatrix, mult, lenght, lenght);

	delete[] mult;

	return reflectionMatrix;
}

/*void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	ofstream file;
	file.open("result.txt");

	file.precision(10);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file << matrix[i * numOfVariables + j] << "     ";
		}
		file << endl;
	}
	file << "--------------------------------------";
	file.close();
}*/

double *GetMatrixForEachStep(double *reflectionMatrix, int numOfColumns, int numOfRows, int step) {
	double *resultMatrix = GetIdentityMatrix(numOfColumns, numOfRows);

	for (int i = step; i < numOfRows; i++) {
		for (int j = step; j < numOfColumns; j++) {
			resultMatrix[i * numOfColumns + j] = reflectionMatrix[(i - step) * (numOfColumns - step) + (j - step)];
		}
	}

	return resultMatrix;
}

double *BackElimination(double *matrix, double *constTermsVec, int numOfVariables, int numOfEquations) {
	double *solutions = new double[numOfVariables];
	int numOfSolutions = numOfVariables - 1;
	int curRow = numOfEquations - 1;
	int curColumn = numOfVariables - 1;

	while (curRow >= 0) {
		double sum = 0;
		for (int i = curColumn + 1; i < numOfVariables; i++) {
			sum += matrix[curRow * numOfVariables + i] * solutions[i];
		}
		solutions[numOfSolutions] = (constTermsVec[curRow] - sum) / matrix[curRow * numOfVariables + curColumn];
		curRow--;
		curColumn--;
		numOfSolutions--;
	}
	
	return solutions;
}

double *ReflectionMethod(double *matrixA, double *matrixB, int numOfVariables, int numOfEquations) {

	//цикл по столбцам
	for (int i = 0; i < numOfVariables - 1; i++) {
		double *normal = Normal(matrixA, numOfVariables, numOfEquations, i);
		double *reflectionMatrix = ReflectionMatrix(normal, numOfEquations - i);
		double *factorMatrix = GetMatrixForEachStep(reflectionMatrix, numOfVariables, numOfEquations, i);
		double *newMatrixA = MultMatrixAB(factorMatrix, numOfVariables, numOfEquations, matrixA, numOfVariables, numOfEquations);
		double *newMatrixB = MultMatrixAB(factorMatrix, numOfVariables, numOfEquations, matrixB, 1, numOfEquations);
		
		delete[] matrixA;
		delete[] matrixB;

		matrixA = newMatrixA;
		matrixB = newMatrixB;

		delete[] normal;
		delete[] reflectionMatrix;
		delete[] factorMatrix;
	}
	double *matrixX = BackElimination(matrixA, matrixB, numOfVariables, numOfEquations);
	delete[] matrixA;
	delete[] matrixB;

	return matrixX;
}

double* Check(double *X, double *A, int numOfVariables, int numOfEqiations) {
	double *res = new double[numOfEqiations];

	for (int i = 0; i < numOfEqiations; i++) {
		res[i] = 0;
	}

	for (int i = 0; i < numOfEqiations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			res[i] += X[j] * A[i * numOfVariables + j];
		}
	}

	return res;
}

/*int main(void) {
	ifstream file;
	file.open("matr1.txt");

	double *dataA;
	int numOfVariablesA, numOfEquationsA;
	
	file >> numOfEquationsA >> numOfVariablesA;

	dataA = new double[numOfVariablesA * numOfEquationsA];

	//считыванием матрицу коэффициентов
	for (int i = 0; i < numOfEquationsA; i++) {
		for (int j = 0; j < numOfVariablesA; j++) {
     		file >> dataA[i * numOfVariablesA + j];
		}
	}
	//--------------------------------------------------------------------------------
	double *dataB;
    dataB = new double[1 * numOfEquationsA];

	//считыванием матрицу коэффициентов
	for (int i = 0; i < numOfEquationsA; i++) {
	  file >> dataB[i];
    }

	double *matrixX = ReflectionMethod(dataA, dataB, numOfVariablesA, numOfEquationsA);
	PrintMatrix(matrixX, 1, numOfEquationsA);

	delete[] matrixX;
	_CrtDumpMemoryLeaks();
}*/



