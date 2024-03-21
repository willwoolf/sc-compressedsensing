#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(initializer_list<double> l) : v(l) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) 
	{ 
		return v[index];
	}

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		return v[index]; 
	}

	// operator overloads
	MVector operator-();

	friend MVector operator*(const double& lhs, const MVector& rhs);
	friend MVector operator*(const MVector& lhs, const double& rhs);
	friend MVector operator+(const MVector& lhs, const MVector& rhs);
	friend MVector operator-(const MVector& lhs, const MVector& rhs);
	friend MVector operator/(const MVector& lhs, const double& rhs);
	friend bool operator==(const MVector& lhs, const MVector& rhs);

	friend ostream& operator<<(ostream& os, const MVector& v);
	
	// get size
	int size() const { return v.size(); } // number of elements

	//vector norms
	double LInfNorm() const;
	double L2Norm() const;
	
	// threshold
	void threshold(int);
	void initialize_normal(int);


private:
	vector<double> v;
};

// dot product
double dot(const MVector& lhs, const MVector& rhs);

// comparison for sorting
bool abscmp(const double&, const double&);

double rand_normal();




#endif
