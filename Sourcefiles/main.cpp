#include "Vector.h"
#include "Matrix.h"
#include "MatrixCalculation.h"
#include <iostream>
#include <fstream>

MatrixCalculation matoperat;

int main()
{
	Vector* x = new Vector(2);
	Vector* b = new Vector(2);
	Matrix* A = new Matrix(2, 2);
	A->setCoeff(0, 0, 2);
	A->setCoeff(0, 1, 3);
	A->setCoeff(1, 0, 4);
	A->setCoeff(1, 1, 5);

	b->setCoeff(0, 1.5);
	b->setCoeff(1, 2);
	b->print();
	matoperat.PLUsolve(A, b, x);
	x->print();
	std::ofstream of("Matrix2.txt");
	A->print(of);

	return 0;
}
