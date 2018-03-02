#include <iostream>
#include <iomanip>
#include <fstream>
#include "linearSystemsSolver.h"

using namespace std;

using dbl = double;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage: \nFirst argument - file name with A and b vectors \nSecond argument - name of method to use" << endl;
        return 0;
    }
    string fileName(argv[1]);
    string methodName(argv[2]);
    ifstream fin(fileName);

    int n;
    fin >> n;

    vector<vector<dbl>> a(n, vector<dbl>(n + 1, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> a[i][j];

    for (int i = 0; i < n; i++)
    {
        fin >> a[i][n];
    }

    linearSystemsSolver solver;
    answer_t response = solver.solve(a, methodName.c_str());

    cout << response.status << endl;
    cout << "answer:\n";
    for (auto i : response.solution)
        cout << setprecision(4) << i << endl;
}
