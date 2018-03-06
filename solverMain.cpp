#include <iostream>
#include <iomanip>
#include <fstream>
#include "linearSystemsSolver.h"

using namespace std;
using namespace linearSystemsSolver;

void solveAll(const matrix_t &a)
{
    std::vector<std::string> allMethods = getAllAvailableMethods();
    for (const auto &name : allMethods)
    {
        answer_t ans = solve(a, name.c_str());
        for (auto i : ans.solution)
            cout << i << " ";
        cout << ans.status << " by method: " << name << endl;
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage: \nFirst argument - file name with A and b vectors \nSecond argument - name of method to use\n";
        return 0;
    }
    string fileName(argv[1]);
    string methodName(argv[2]);
    ifstream fin(fileName);

    int n;
    fin >> n;

    using dbl = double;
    matrix_t a(n, vector<dbl>(n + 1, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n; i++)
    {
        fin >> a[i][n];
    }

    if (methodName == "all")
    {
        solveAll(a);
        return 0;
    }

    answer_t response = solve(a, methodName.c_str());

    cout << "Cond of system: " << cond(a) << endl;
    cout << response.status << endl;

    if (response.amountOfSolutions != Inf && response.amountOfSolutions != 0)
    {
        cout << "answer:\n";
        for (auto i : response.solution)
            cout << setprecision(4) << i << endl;
    }
}


