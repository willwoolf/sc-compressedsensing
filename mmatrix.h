#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include <iostream>
#include "mvector.h"

using namespace std;

// Class that represents a mathematical matrix
class MMatrix
{
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	friend MVector operator*(const MMatrix& A, const MVector& v);
	friend MMatrix operator*(const MMatrix& A, const MMatrix& B);
	MMatrix transpose() const;

	void initialize_normal();

	friend ostream& operator<<(ostream& os, const MMatrix& A);

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

private:
	unsigned int nRows, nCols;
	vector<double> A;
};


#endif
