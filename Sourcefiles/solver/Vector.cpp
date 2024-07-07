/*
define vector
*/
#include "Vector.h"

Vector::Vector()
{
	cdp_vecCoeff = nullptr;
}

Vector::Vector(int nr) 
{
	ci_vecRow = nr;
	cdp_vecCoeff = new double[nr];
}

Vector::~Vector() 
{
	if (cdp_vecCoeff) 
	{
		delete[]cdp_vecCoeff;
		cdp_vecCoeff = nullptr;
	}	
}

void Vector::setNumRows(int NumRows) 
{
	ci_vecRow = NumRows;
}

void Vector::setCoeff(double* vecCoeff)
{
	cdp_vecCoeff = vecCoeff;
}

void Vector::zero()
{
	for (int i = 0; i < ci_vecRow; i++)
	{
		cdp_vecCoeff[i] = 0;
	}
}

int Vector::i_getNumRows() const
{
	return ci_vecRow;
}

void Vector::setCoeff(int i, double value)
{
	cdp_vecCoeff[i] = value;
}

void Vector::addCoeff(int i, double value)
{
	cdp_vecCoeff[i] = cdp_vecCoeff[i] + value;
}

double Vector::d_getCoeff(int i)
{
	return cdp_vecCoeff[i];
}

double Vector::get_modulus()const
{
	double sum = 0;
	for (int i = 0; i < ci_vecRow; i++)
	{
		sum = sum + cdp_vecCoeff[i] * cdp_vecCoeff[i];
	}
	return sqrt(sum);
}

void Vector::print(ofstream& fout)
{
	fout << "向量行数 : " << ci_vecRow << std::endl;
	for (int i = 0; i < ci_vecRow; i++)
	{
		fout << setiosflags(ios::scientific) << setprecision(12) << setw(15) << cdp_vecCoeff[i] << std::endl;
	}
}

void Vector::print()
{
	std::cout << "向量行数 : " << ci_vecRow << std::endl;
	for (int i = 0; i < ci_vecRow; i++)
	{
		std:: cout << setiosflags(ios::scientific) << setprecision(12) << setw(15) << cdp_vecCoeff[i] << std::endl;
	}
}


