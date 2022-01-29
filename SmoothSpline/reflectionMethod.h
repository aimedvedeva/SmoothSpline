#pragma once

double *ReflectionMethod(double *matrixA, double *matrixB, int numOfVariables, int numOfEquations);
//double *TransposeMatrix(double *matrix, int numOfColumns, int numOfRows);
double *MultMatrixAB(double *matrixA, int numOfVariablesA, int numOfEquationsA, double *matrixB, int numOfVariablesB, int numOfEquationsB);
