#include"Matrix.h"

Matrix::Matrix(int nr, int nc)
{
	ci_matRow = nr;
	ci_matCol = nc;
	cdp_matCoeff = new double[nr * nc];
};

Matrix::~Matrix()
{
	delete[]cdp_matCoeff;
	cdp_matCoeff = nullptr;
}

void Matrix::zero()
{
	for (int i = 0; i < ci_matRow; i++)
	{
		for (int j = 0; j < ci_matCol; j++)
		{
			cdp_matCoeff[i * ci_matCol + j] = 0;
		}
	}
}

int Matrix::i_getNumRows()const
{
	return ci_matRow;
}

int Matrix::i_getNumCols()const
{
	return ci_matCol;
}

void Matrix::setCoeff(int i, int j, double value)
{
	cdp_matCoeff[i * ci_matCol + j] = value;
}

void Matrix::addCoeff(int i, int j, double value)
{
	cdp_matCoeff[i * ci_matCol + j]=cdp_matCoeff[i * ci_matCol + j] + value;
}

double Matrix::d_getCoeff(int i, int j)
{
	return cdp_matCoeff[i * ci_matCol + j];
}

void  Matrix::print(ofstream& fout) const
{
	fout << setiosflags(ios::scientific)
		<< setprecision(12);
	fout << "Matrix Size: " << "(" << ci_matRow << ", " << ci_matCol << ")" << endl;
	for (int i = 0; i < ci_matRow; i++)
	{
		//fout << "Row(" << i << "):" << endl;
		for (int j = 0; j < ci_matCol; j++)
		{
			fout << setw(15) << cdp_matCoeff[i * ci_matCol + j] << '\t';
		}
		fout << endl;
	}
}