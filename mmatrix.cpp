#include <vector>
#include <iostream>
#include <iomanip>
#include <ios>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "mmatrix.h"
#include "mvector.h"

using namespace std;

MVector operator*(const MMatrix& A, const MVector& v)
{
    // check compatibility
    if (A.Cols() != v.size())
    {
        cout << "wrong dims" << endl;
        exception e;
        throw e;
    }
    // initialise output vector
    int vout_size = A.Rows();
    MVector vout = MVector(vout_size);

    for (int i = 0; i < A.Rows(); i++)
    {
        // sum across j
        double sum = 0;
        for (int j = 0; j < A.Cols(); j++)
        {
            sum += A(i, j)*v[j];
        }
        vout[i] = sum;        
    }
    return vout;
}

MMatrix operator*(const MMatrix& A, const MMatrix& B)
{
    // compatibility check
    assert(A.nCols == B.nRows);

    // construct the return matrix of correct size
    MMatrix C(A.nRows, B.nCols);

    // compute all C(i, j)
    for (int i = 0; i < C.nRows; i++)
    {
        for (int j = 0; j < C.nCols; j++)
        {
            double sum = 0;
            for (int k = 0; k < A.nCols; k++)
            {
                sum += A(i, k)*B(k, j);
            }
            C(i, j) = sum;
        }
    }
    return C;
}

MMatrix MMatrix::transpose() const
{
    MMatrix C(nCols, nRows);
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            // hacky. i and j are rows and colums respectively of A, opposite for C
            C(j, i) = A[j + i*nCols];
        }
    }
    return C;
}

void MMatrix::initialize_normal()
{
    for (vector<double>::iterator term = A.begin(); term != A.end(); term++)
    {
        // set every entry to output of a zero mean, variance 1/m gaussian
        *term = rand_normal()/sqrt(nCols);
    }
    
}

ostream& operator<<(ostream& os, const MMatrix& A)
{
    for (int i = 0; i < A.Rows(); i++)
    {
        os << "|  ";
        for (int j = 0; j < A.Cols(); j++)
        {
            os << setw(3) << A(i, j) << "  ";
        }
        os << "|\n";
    }
    return os;
}