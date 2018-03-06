#include <vector>
#include <cmath>
#include "linearSystemsSolver.h"

typedef double dbl;

int gauss(std::vector < std::vector<dbl> > a, std::vector<dbl> &ans)
{
    int amountOfEq = (int) a.size();
    int amountOfVars = (int) a[0].size() - 1;

    std::vector<int> placeHolder(amountOfVars, -1);
    for (int col = 0, row = 0; (col < amountOfVars) && (row < amountOfEq); ++col)
    {
        int sel = row;

        for (int i = row; i < amountOfEq; ++i)
        {
            if (std::fabs(a[i][col]) > std::fabs(a[sel][col]))
                sel = i;
        }
        if (std::fabs(a[sel][col]) < linearSystemsSolver::EPS) continue;
        for (int i = col; i <= amountOfVars; ++i) std::swap(a[sel][i], a[row][i]);

        placeHolder[col] = row;

        for (int i = 0; i < amountOfEq; ++i)
        {
            if (i != row)
            {
                dbl c = a[i][col] / a[row][col];
                for (int j = col; j <= amountOfVars; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        }
        ++row;
    }

    ans.assign(amountOfVars, 0);
    for (int i = 0; i < amountOfVars; ++i)
        if (placeHolder[i] != -1)
            ans[i] = a[placeHolder[i]][amountOfVars] / a[placeHolder[i]][i];
    for (int i = 0; i < amountOfEq; ++i)
    {
        dbl sum = 0;
        for (int j = 0; j < amountOfVars; ++j) sum += ans[j] * a[i][j];
        if (std::fabs(sum - a[i][amountOfVars]) > linearSystemsSolver::EPS) return 0;
    }
    for (int i = 0; i < amountOfVars; ++i) if (placeHolder[i] == -1) return linearSystemsSolver::Inf;
    return 1;
}
