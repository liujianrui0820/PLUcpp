#pragma once

#include"Matrix.h"
#include<iomanip>
#include<iostream>
#include<fstream>
using namespace std;

class Matrix
{
private:
	double* cdp_matCoeff;	//coefficient of matrix
	int ci_matRow, ci_matCol;	//The row and column positions of a matrix
	Matrix();	//Prevent defining a matrix with uncertain rows and columns
	//Matrix(const Matrix& o_matrix);	Prevent matrix objects from being copied

public:
	//constructor and destructor
	Matrix(int nr, int nc);
	~Matrix();

	//function
	void zero();	//Initialize Zero Matrix
	int i_getNumRows()const;	//Obtain the number of matrix rows
	int i_getNumCols()const;	//Obtain the number of matrix cols
	void setCoeff(int i, int j, double value);	//set values to i-th row and j-th column
	void addCoeff(int i, int j, double value);	//add values to i-th row and j-th column
	double d_getCoeff(int i, int j);	//get cdp_matCoeff[i][j];

	void print(ofstream& fout) const;

};