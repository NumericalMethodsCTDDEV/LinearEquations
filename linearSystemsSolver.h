#ifndef LINEARSYSTEMSSOLVER_H
#define LINEARSYSTEMSSOLVER_H

#include <string>
#include <vector>
#include <limits.h>

using matrix_t = std::vector<std::vector<double>>;

namespace linearSystemsSolver
{
    const double EPS =  0.00000001;

    struct answer_t
    {
        std::string status;
        std::vector<double> solution;
        int amountOfSolutions;
    };

    const int Inf = INT_MAX - 100;

    answer_t solve(const matrix_t &system, const char *methodName);

    std::vector<std::string> getAllAvailableMethods();

    matrix_t transpose(const matrix_t &);

    double cond(const matrix_t &);

    double determinant(const matrix_t &);

    double norm(const matrix_t&);
}

#endif // LINEARSYSTEMSSOLVER_H
