/*
Vector class

*/
#pragma once

#include<fstream>
#include<cmath>
#include<iomanip>
#include<iostream>

using namespace std;

class Vector
{
public:
	double* cdp_vecCoeff;	//Class Double Point vector Coeffient
	int ci_vecRow;	// class int row of vector

public:
	//	Constructor and Destructor	构造函数和析构函数
	Vector();	//construct a vector(default)
	Vector(int nr);	//	construct a vector (Number of rows: n)
	~Vector();		//	release dynamic memory

	//function
	void setNumRows(int NumRows);
	void setCoeff(double* vecCoeff);

	void zero();	//initialiaze zero vector
	int i_getNumRows() const;	//get NumRows
	void setCoeff(int i, double value);	//set value of cdp_vecCoeff[i]
	void addCoeff(int i, double value);	//add value to cdp_vecCoedd[i]
	double d_getCoeff(int i);	//get values of cdp_vecCoeff[i]
	double get_modulus()const;//get modulus of vector
	void print(ofstream& fout);	//fout
	void print();	//cout

private:
	// never to be used copy constructor
	Vector(const Vector& o_vector);
};