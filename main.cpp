#include <iostream>
#include <fstream>
#include <thread>

#include "mvector.h"
#include "mmatrix.h"

using namespace std;

// non-class functions
int SDLS(const MMatrix& A, const MVector& b, MVector& x, int maxIterations, double tol, bool writesFile = false, string filename = "")
{
    // r is the residual for solving the normal equations
    MVector r = A.transpose() * (b - A*x);

    // block for writing data. f must be created but if !writesFile we write nothing
    ofstream f;
    if (writesFile && filename != "")
    {
        f.open(filename);
        if (!f)
        {
            exception e;
            throw e;
        }
        f << "Steepest Descent Algorithm implementation" << "\n";
        f << "A = " << "\n" << A << "\n";
        f << "b = " << b << "\n";
        f << "Iterations for x:\n";
        f << x << "\n";
    }

    int iter = 0;
    while (iter < maxIterations && r.L2Norm() > tol)
    {
        // perform the iteration 
        double alpha = (dot(r, r))/(dot(A*r, A*r));
        x = x + alpha * r;
        r = r - alpha*(A.transpose()*(A*r));

        if (writesFile && filename != "")
        {
            f << x << "\n";
        }
        
        iter++;
    }
    if (writesFile && filename != "")
    {
        f.close();
    }
    
    if (r.L2Norm() > tol)
    {
        // iteration fails to converge to a solution
        return 0;
    }

    return iter;
}

int NIHT(const MMatrix& A, const MVector& b, MVector&x, int k, const int& maxIterations, const double& tol, const bool& writesFile = false, string filename = "")
{
    // compute vector norms and sqrt(tol) for convergence checking
    double x_infnorm = x.LInfNorm(), x_2norm = x.L2Norm();
    double laxtol = sqrt(tol);
    
    // Initialise starting vector
    x = A.transpose()*b;
    x.threshold(k); // Get k largest values
    MVector r = A.transpose() * (b - A*x);

    // block for writing data
    ofstream f;
    if (writesFile && filename != "")
    {
        f.open(filename);
        if (!f)
        {
            exception e;
            throw e;
        }
        f << "Normalised Iterative Hard Thresholding Algorithm implementation" << "\n";
        f << "A is a random " << A.Rows() << " by " << A.Cols() << " matrix.\n";
        f << "x is a random " << x.size() << "-vector with sparsity " << k << ".\n";
        f << "Norm of the residual:\n";
        f << r.L2Norm() << "\n";

        // f << "Iterations for x:\n";
        // f << x << "\n";
    }
    
    // begin iteration
    int iter = 0;
    while (iter < maxIterations && r.L2Norm() > tol)
    {
        double alpha = (dot(r, r))/(dot(A*r, A*r));
        // compute and threshold x
        MVector xPrev = x;
        x = x + alpha*r;
        x.threshold(k);

        // cout << x << "\n";

        // compute residual
        r = A.transpose() * (b - A*x);    

        // do some convergence testing
        if (abs(xPrev.LInfNorm() - x.LInfNorm()) < tol && r.L2Norm() > laxtol)
        {
            return 0;
        }
        
        if (writesFile && filename != "")
        {
            f << r.L2Norm() << "\n";
            // f << x << "\n";
        }
        iter++;
    }

    if (writesFile && filename != "")
    {
        f.close();
    }
    if (r.L2Norm() > tol)
    {
        // iteration has failed to converge
        return 0;
    }
    return iter;

}

double stability(int m, int n, int k, int T)
{
    /*
    Stability returns the probability p(m) of successful recovery of NIHT given dimensions, sparsity, and the number of times T to perform NIHT.
    */
    double successes = 0;
    for (int i = 0; i < T; i++)
    {
        // make a k-sparse vector x, make A, form b=Ax
        MVector x(n);
        x.initialize_normal(k);
        MMatrix A(m, n);
        A.initialize_normal();
        MVector b = A*x;

        double x_linfnorm_orig = x.LInfNorm();

        int iterations = NIHT(A, b, x, k, 100000, 1e-6);
        if (abs(x.LInfNorm()-x_linfnorm_orig) < 1e-3 && iterations > 0)
        {
            successes++;
        }
    }
    return successes/double(T);
}

int main()
{
    srand(time(NULL));
    
    // (1) generating a matrix and testing that transpose() works
    // MMatrix A(2, 3);
    // A(0, 0) = 1;
    // A(1, 1) = 2;
    // A(2, 0) = 1;

    // cout << A << endl;
    // A.transpose();
    // cout << A << endl;

    // (2) least squares convergence
    // MMatrix A(3, 2);
    // A(0, 0) = 1;
    // A(1, 0) = 2;
    // A(2, 0) = -1;
    // A(0, 1) = 2;
    // A(1, 1) = 1;
    // A(2, 1) = 0;

    // MVector b = {10, -1, 0};
    
    // MVector x = {0, 0};
    // cout << A << endl;
    // int n = SDLS(A, b, x, 1000, 1e-6, true, "sdsolveA1.txt");

    // A(2, 0) = 1.8;
    // A(2, 1) = -2;
    // x = {0, 0};
    // n = SDLS(A, b, x, 1000, 1e-6, true, "sdsolveA2.txt");

    // A(2, 0) = -2;
    // A(2, 1) = -2;
    // x = {0, 0};
    // n = SDLS(A, b, x, 1000, 1e-6, true, "sdsolveA3.txt");
    // cout << n << ", " << x << ", " << A*x << endl;

    // (3) For generating a histogram of the entries for a 1000x1 random matrix
    // MMatrix N(1000, 1);
    // N.initialize_normal();

    // ofstream f;
    // f.open("normaldisttest.txt");
    // if (!f)
    // {
    //     exception e;
    //     throw e;
    // }
    // f << N << endl;
    // f.close();
    
    // (4) perform normalised iterative hard thresholding for a single example
    double eqnratio, sparsity;
    eqnratio = 0.6, sparsity = 0.4;
    int n = 5;
    int m = eqnratio * n;
    int k = sparsity * n;
    MMatrix A(m, n);
    MVector b(m);
    MVector x(n);

    A.initialize_normal();
    x.initialize_normal(k);

    b = A*x;

    cout << A << endl;
    cout << x << endl;
    cout << b << endl;

    int iterations = NIHT(A, b, x, k, 100000, 1e-6);
    cout << "iteration NIHT performed in " << iterations << " steps\n";
    cout << x << endl;

    // (5) compute the phase transition
    // int n = 200, k = 20, T = 50;

    // for (int m = 4; m <= 199; m += 5)
    // {
    //     double p = stability(m, n, k, T);
    //     cout << "For m = " << m << ", n = " << n
    //     << ", the probability of recovery for sparsity " << k << " is " << p << endl;
    // }

    // (6) testing
    // cout << stability(60, 50, 10, 50) << "\n";
    // cout << stability(120, 100, 20, 50) << "\n";
    // cout << stability(240, 200, 40, 50) << "\n";

    // MMatrix A(3,2);
    // int k = 2;
    // A.initialize_normal();
    // MVector xstar(2);
    // xstar.initialize_normal(2);
    // MVector xiter(2), yiter(2);

    // cout << SDLS(A, A*xstar, xiter, 1000, 1e-4, true, "randomsdlstest.txt") << endl;
    // cout << NIHT(A, A*xstar, yiter, k, 1000, 1e-4, true, "randomnihttest.txt") << endl;

    // (7) 2d figure generator
    /*
    compute the figure as a matrix.
    n = 200,
    m ranging 4:5:199 is 5*(1:1:40)-1
    k ranging 5:5:100 is 5*(1:1:20)
    */

    // ofstream f;
    // string filename = "placeholder.txt";
    // f.open(filename);

    // // generate a results matrix with all the values and output this to a file.
    // // generate a file for each value of m with the results

    // int n = 200;
    // int T = 50;
    // int rows = 40, cols = 20;
    // for (int m = 1; m <= rows; m++)
    // {
    //     int M = 5 * m - 1;
    //     for (int k = 1; k <= cols; k++)
    //     {
    //         int K = 5*k;
    //         double p = stability(M, n, K, T);
            
    //         cout << "For m = " << M << ", n = " << n
    //         << ", the probability of recovery for sparsity " << K << " is " << p << endl;
            
    //         // index A by row m, column k
    //         f << "A(" << M << ", " << K << ") = " << p << "\n";
        
    //     }
    // }
    // f.close();
    
    return 0;
}