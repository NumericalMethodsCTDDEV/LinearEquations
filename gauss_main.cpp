#include <iostream>
#include <iomanip>
#include "gauss.h"
#include <fstream>

using namespace std;

using dbl = double;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: pass file name with A and b vectors to arguments" << endl;
    }
    string fileName(argv[1]);
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

    vector<dbl> ans;
    int status = gauss(a, ans);

    cout << "amount of solutoins: " << status << endl;
    cout << "answer:\n";
    for (auto i : ans)
        cout << setprecision(4) << i << endl;
}
