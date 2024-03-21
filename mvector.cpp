#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <cassert>
#include <ctime>

#include "mvector.h"

using namespace std;

// constructors in header

MVector operator*(const double& lhs, const MVector& rhs)
{
    MVector temp = rhs;
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] *= lhs;
    }
    return temp;
}

MVector operator*(const MVector& lhs, const double& rhs)
{
    return rhs*lhs;
}

MVector operator+(const MVector& lhs, const MVector& rhs)
{
    if (lhs.size() != rhs.size())
    {
        exception e;
        throw e;
    }
    MVector temp(lhs.size());
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] = lhs[i] + rhs[i];
    }
    return temp;
}

MVector operator-(const MVector& lhs, const MVector& rhs)
{
    return lhs+(-1*rhs);
}

// alternative overload such that given MVector v, -v is allowed
MVector MVector::operator-()
{
    MVector vec(v.size());
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] = -v[i];
    }
    return vec;
}

MVector operator/(const MVector& lhs, const double& rhs)
{
    return (1/rhs)*lhs;
}

bool operator==(const MVector& lhs, const MVector& rhs)
{
    assert(lhs.size() == rhs.size());
    for (int i = 0; i < lhs.size(); i++)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }
    }
    return true;
    
}

ostream& operator<<(ostream& os, const MVector& v)
{
	int n = v.size();
	os << "(";
	for (int i = 0; i < n-1; i++)
    {
        os << setw(10) << v[i] << ", ";
    }
    os << v[n-1];
	os << ")";
	return os;
}

double MVector::LInfNorm() const
{
    double maxAbs = 0;
    size_t s = size();
    for (int i=0; i<s; i++)
    {
        maxAbs = max(abs(v[i]), maxAbs);
    }
    return maxAbs;
}

double MVector::L2Norm() const
{
    double sum = 0;
    for (int i = 0; i < v.size(); i++)
    {
        sum += v[i]*v[i];
    }
    return sqrt(sum);
}

double dot(const MVector& lhs, const MVector& rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw invalid_argument("error: dot product not defined for vectors of non-equal dimension");
    }
    double sum = 0;
    for (int i = 0; i < lhs.size(); i++)
    {
        sum += lhs[i]*rhs[i];
    }
    return sum;
}

bool abscmp(const double& a, const double& b)
{
    return abs(a) < abs(b);
}

void MVector::threshold(int k)
{
    vector<double> v_copy = v;
    // sort a copy by absolute value, then identify the k-th largest value
    sort(v_copy.begin(), v_copy.end(), abscmp);
    // this picks out the k-th largest value
    double bound = v_copy[v_copy.size() - k];

    // iterate through the MVector contents performing the thresholding
    for (vector<double>::iterator term = v.begin(); term != v.end(); term++)
    {
        if (abscmp(*term, bound))
        {
            *term = 0;
        }   
    }
}

void MVector::initialize_normal(int k)
{
    for (vector<double>::iterator term = v.begin(); term != v.end(); term++)
    {
        *term = rand_normal();
    }

    // complete copy of threshold
    vector<double> v_copy = v;
    sort(v_copy.begin(), v_copy.end(), abscmp);
    double bound = v_copy[v_copy.size() - k];
    for (vector<double>::iterator term = v.begin(); term != v.end(); term++)
    {
        if (abscmp(*term, bound))
        {
            *term = 0;
        }   
    }   
}

double rand_normal()
{
    static const double pi = 3.141592653589793238;
	double u = 0;
	while (u == 0) // loop to ensure u nonzero, for log
	{
        u = rand() / static_cast<double>(RAND_MAX);
    }
    double v = rand() / static_cast<double>(RAND_MAX);
	return sqrt(-2.0*log(u))*cos(2.0*pi*v);
}