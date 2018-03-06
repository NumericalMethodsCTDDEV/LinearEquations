#include "linearSystemsSolver.h"
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>

using vector_t = std::vector<double>;

int gauss(matrix_t, vector_t &);
int jacobi(matrix_t, vector_t &);
int seidel(matrix_t, vector_t &);
int sequantialRelaxation(matrix_t, vector_t &);
int descent(matrix_t, vector_t &);

using solveMethod_t = int (matrix_t, vector_t &);

static std::unordered_map<std::string, solveMethod_t *> nameToMethod =
{
    {"gauss", &gauss},
    {"jacobi", &jacobi},
    {"seidel", &seidel},
    {"sequantialRelaxation", &sequantialRelaxation},
    {"descent", &descent}
};

static std::string parseStatus(int st)
{
    switch (st)
    {
    case 0:
        return "No solutions";
    case 1:
        return "One solution";
    case linearSystemsSolver::Inf:
        return "Inf solutions";
    default:
        return "Something went wrong";
    }
}

static inline void getCofactor(const matrix_t &mat, matrix_t &cofactors, int p, int q, int n)
{
    int i = 0, j = 0;
    for (int row = 0; row < n; row++)
        for (int col = 0; col < n; col++)
            if (row != p && col != q)
            {
                cofactors[i][j++] = mat[row][col];
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
}

static double getDeterminant(const matrix_t &mat, int n)
{
    double d = 0;
    if (n == 1)
        return mat[0][0];
    matrix_t cofactors(mat.size(), std::vector<double>(mat.size(), 0));
    int sign = 1;
    for (int f = 0; f < n; f++)
    {
        getCofactor(mat, cofactors, 0, f, n);
        d += sign * mat[0][f] * getDeterminant(cofactors, n - 1);
        sign = -sign;
    }
    return d;
}

static double _norm(const matrix_t &a, bool needMax = true)
{
    size_t n = a.size();
    n = std::min(n, a[0].size());
    double sMax = (needMax ? -1 : std::numeric_limits<double>::max());
    for (size_t i = 0; i < n; ++i)
    {
        double s = 0;
        for (size_t j = 0; j < n; ++j)
            s += std::fabs(a[j][i]);
        if (needMax) sMax = std::max(sMax, s);
        else sMax = std::min(sMax, s);
    }
    return sMax;
}

namespace linearSystemsSolver
{
    answer_t solve(const matrix_t &system, const char *methodName)
    {
        const std::string name(methodName);
        vector_t ans;
        if (!nameToMethod.count(name))
            return {"No such method", {}, 0};
        const solveMethod_t *f = nameToMethod[name];
        int amount = f(system, ans);
        return {parseStatus(amount), ans, amount};
    }

    matrix_t transpose(const matrix_t &system)
    {
        size_t n = system.size();
        n = std::min(n, system[0].size());
        matrix_t s = system;
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i; j < n; ++j)
                std::swap(s[i][j], s[j][i]);
        return s;
    }

    double determinant(const matrix_t &matr)
    {
        return getDeterminant(matr, static_cast<int>(matr.size()));
    }

    double cond(const matrix_t &a)
    {
        return _norm(a) / _norm(a, false);
    }

    double norm(const matrix_t &a)
    {
        return _norm(a);
    }

}
