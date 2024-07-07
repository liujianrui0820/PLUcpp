#include"MatrixCalculation.h"
#include<iostream>
void MatrixCalculation::matTranspose(Matrix* A, Matrix* A_trans)
{
	int nr = A->i_getNumRows();
	int nc = A->i_getNumCols();
	double temp;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			temp = A->d_getCoeff(i, j);
			A_trans->setCoeff(j, i, temp);
		}
	}
}

void MatrixCalculation::matMultiply(Matrix* A, Matrix* B, Matrix* C)
{
	int A_nr = A->i_getNumRows();
	int A_nl = A->i_getNumCols();
	int B_nr = B->i_getNumRows();
	int B_nl = B->i_getNumCols();

	for (int i = 0; i < A_nr; i++)
	{
		double temp = 0;
		for (int j = 0; j < B_nl; j++)
		{
			double temp = 0;
			for (int k = 0; k < B_nr; k++)
			{
				temp = (A->d_getCoeff(i, k)) * (B->d_getCoeff(k, j)) + temp;
			}
			C->setCoeff(i, j, temp);
		}
	}
}


void MatrixCalculation::matMultiply(Matrix* A, Vector* b, Vector* c)
{
	int A_nr = A->i_getNumRows();
	int A_nl = A->i_getNumCols();
	for (int i = 0; i < A_nr; i++)
	{
		double temp = 0;
		for (int j = 0; j < A_nl; j++)
		{
			temp = (A->d_getCoeff(i, j)) * (b->d_getCoeff(j)) + temp;
		}
		c->setCoeff(i, temp);
	}
}

void MatrixCalculation::matMultiply(Vector* a, Vector* b, double& c)
{
	int nr = a->i_getNumRows();
	c = 0;
	for (int i = 0; i < nr; i++)
	{
		c = (a->d_getCoeff(i)) * (b->d_getCoeff(i)) + c;
	}
}

void MatrixCalculation::matMultiply(Matrix* A, double b, Matrix* C)
{
	int A_nr = A->i_getNumRows();
	int A_nl = A->i_getNumCols();
	double temp;
	for (int i = 0; i < A_nr; i++)
	{
		for (int j = 0; j < A_nl; j++)
		{
			temp = A->d_getCoeff(i, j) * b;
			C->setCoeff(i, j, temp);
		}
	}
}

void MatrixCalculation::matMultiply(Vector* a, double b, Vector* c)
{
	int nr = a->i_getNumRows();
	double temp;
	for (int i = 0; i < nr; i++)
	{
		temp = a->d_getCoeff(i) * b;
		c->setCoeff(i, temp);
	}
}

void MatrixCalculation::matAdd(Matrix* A, Matrix* B, Matrix* C)
{
	int A_nr = A->i_getNumRows();
	int A_nc = A->i_getNumCols();
	int B_nr = B->i_getNumRows();
	int B_nc = B->i_getNumCols();
	if (A_nr != B_nr || A_nc != B_nc)
	{
		throw std::invalid_argument("Matrix dimension mismatch");
		exit(0);
	}
	for (int i = 0; i < A_nr; i++)
	{
		for (int j = 0; j < A_nc; j++)
		{
			double temp = A->d_getCoeff(i, j) + B->d_getCoeff(i, j);
			C->setCoeff(i, j, temp);
		}
	}
}

void MatrixCalculation::matAdd(Vector* a, Vector* b, Vector* c)
{
	int a_nr = a->i_getNumRows();
	int b_nr = b->i_getNumRows();
	if (a_nr != b_nr)
	{
		throw std::invalid_argument("Vector dimension mismatch");
		exit(0);
	}
	for (int i = 0; i < a_nr; i++)
	{
		double temp = a->d_getCoeff(i) + b->d_getCoeff(i);
		c->setCoeff(i, temp);
	}
}

void MatrixCalculation::matMinus(Matrix* A, Matrix* B, Matrix* C)
{
	int A_nr = A->i_getNumRows();
	int A_nc = A->i_getNumCols();
	int B_nr = B->i_getNumRows();
	int B_nc = B->i_getNumCols();
	if (A_nr != B_nr || A_nc != B_nc)
	{
		throw std::invalid_argument("Matrix dimension mismatch");
		exit(0);
	}
	for (int i = 0; i < A_nr; i++)
	{
		for (int j = 0; j < A_nc; j++)
		{
			double temp = A->d_getCoeff(i, j) - B->d_getCoeff(i, j);
			C->setCoeff(i, j, temp);
		}
	}
}

void MatrixCalculation::matMinus(Vector* a, Vector* b, Vector* c)
{
	int a_nr = a->i_getNumRows();
	int b_nr = b->i_getNumRows();
	if (a_nr != b_nr)
	{
		throw std::invalid_argument("Vector dimension mismatch");
		exit(0);
	}
	for (int i = 0; i < a_nr; i++)
	{
		double temp = a->d_getCoeff(i) - b->d_getCoeff(i);
		c->setCoeff(i, temp);
	}
}

//PLU decomposition

void MatrixCalculation::PLUsolve(Matrix* A, Vector* b, Vector* x)
{
	double temp, temp1;
	int Dim;
	int A_r = A->i_getNumRows();
	int A_l = A->i_getNumCols();
	if (A_r != A_l)
	{
		std::cout << "error" << std::endl;
		throw std::invalid_argument("Vector dimension mismatch");
		exit(0);
	}
	else
	{
		Dim = A_r;
	}
	int* P;
	P = new int[Dim];
	Matrix* L, * U;
	L = new Matrix(Dim, Dim);
	U = new Matrix(Dim, Dim);

	PLUdecomposition(A, L, U, P);
	//LUx=b		1,Ly=b;		2，Ux=y
	double* y;
	y = new double[Dim];
	for (int i = 0; i < Dim; i++)
	{
		y[i] = b->d_getCoeff(P[i]);//y[i]=b,	(matrix P),		y[i]=b-coefficient*y
		for (int j = 0; j < i; j++)
		{
			y[i] = y[i] - L->d_getCoeff(i, j) * y[j];
		}
	}
	//	Ux=y
	for (int i = Dim -1; i >=0; i--)
	{
		x->setCoeff(i, y[i]);
		for (int j = Dim - 1; j > i; j--)
		{
			temp = x->d_getCoeff(i) - U->d_getCoeff(i, j) * x->d_getCoeff(j);
			x->setCoeff(i, temp);
		}
		temp1 = x->d_getCoeff(i) / U->d_getCoeff(i, i);
		x->setCoeff(i, temp1);
	}

	delete L;
	delete U;
	delete[] y;
	delete[] P;
	L = nullptr;
	U = nullptr;
	P = nullptr;
	y = nullptr;

}


void MatrixCalculation::PLUdecomposition(Matrix* A, Matrix* L, Matrix* U, int* P)
{
	L->zero();
	U->zero();
	int matDim = A->i_getNumRows();
	int row = 0;
	int temp1;
	double temp2, temp3;

	//Each row is given an initial number to describe the P-matrix
	for (int i = 0; i < matDim; i++)
	{
		P[i] = i;
	}

	double u, l;
	for (int i = 0; i < matDim; i++)
	{
		double PartialPivot = 0;
		for (int j = i; j < matDim; j++)
		{
			if (fabs(A->d_getCoeff(j, i)) > PartialPivot)
			{
				PartialPivot = A->d_getCoeff(j, i);
				row = j;
			}
		}

		if (PartialPivot == 0)
		{
			std::cout << "Matrix Singularity" << std::endl;
			return;
		}
		//Exchange the i-th and j-th rows of the matrix
		temp1 = P[i];
		P[i] = P[row];
		P[row] = temp1;

		temp2 = 0;
		for (int k = 0; k < matDim; k++)
		{
			//把第i行第k列的数存到临时，给第i行第k列赋第row行第j列的值，给第row行第j列赋临时值
			temp2 = A->d_getCoeff(i, k);
			A->setCoeff(i, k, A->d_getCoeff(row, k));
			A->setCoeff(row, k, temp2);
		}

		u = A->d_getCoeff(i, i);
		l = 0;
		for (int j = i+1; j < matDim; j++)		//Traverse each row
		{
			l = A->d_getCoeff(j, i) / u;
			A->setCoeff(j, i, l);
			for (int k = i+1; k < matDim; k++)	//Traverse each column
			{
				temp3 = A->d_getCoeff(j, k) - A->d_getCoeff(i, k) * l;
				A->setCoeff(j, k, temp3);
			}
		}

		//Calculate the L and U matrices
	   //construct U and L
		for (int i = 0; i < matDim; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if (i != j)
				{
					L->setCoeff(i, j, A->d_getCoeff(i, j));
				}
				else
				{
					L->setCoeff(i, j, 1.0);
				}
			}
			for (int k = i; k < matDim; k++)
			{
				U->setCoeff(i, k, A->d_getCoeff(i, k));
			}
		}


	}
}
