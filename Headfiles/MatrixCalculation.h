#pragma once

#include"Vector.h"
#include"Matrix.h"
#include<vector>
#include<cmath>

class MatrixCalculation
{
public:
	//Calculate the transpose of matrix A, and the result is represented by A_trans
	void matTranspose(Matrix* A, Matrix* A_trans);
	//Multiplying matrix A by matrix B, the result is represented by C
	void matMultiply(Matrix* A, Matrix* B, Matrix* C);
	//Multiplying matrix A by vector b, the result is represented by c
	void matMultiply(Matrix* A, Vector* b, Vector* c);
	//Multiplying vector a by vector b, the result is represented by c
	void matMultiply(Vector* a, Vector* b, double& c);
	//Multiplying matrices and vectors by constants
	void matMultiply(Matrix* A, double b, Matrix* C);
	void matMultiply(Vector* a, double b, Vector* c);

	void matAdd(Matrix* A, Matrix* B, Matrix* C);
	void matAdd(Vector* a, Vector* b, Vector* c);
	void matMinus(Matrix* A, Matrix* B, Matrix* C);
	void matMinus(Vector* a, Vector* b, Vector* c);

	//Using PLU decomposition to solve linear systems of equations Ax=b
	void PLUsolve(Matrix* A, Vector* b, Vector* x);

	//void matAdd(Vec)

private:
	void PLUdecomposition(Matrix* A, Matrix* L, Matrix* U, int* P);

};
