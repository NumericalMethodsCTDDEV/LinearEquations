#include <bits/stdc++.h>
#include "linearSystemsSolver.h"
#include <functional>

using namespace std;

const double N = 100;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

void generate(size_t n, std::function<double (int, int, double)> &&genNumber)
{
    std::uniform_real_distribution<> dis(1, N);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            cout << genNumber(i, j, dis(gen)) << " ";
        }
        cout << endl;
    }
    for (size_t i = 0; i < n; ++i)
        cout << dis(gen) << endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3) return 0;

    size_t n = strtoull(argv[1], nullptr, 10);
    cout << n << endl;
    string target(argv[2]);

    if (target == "good")
        generate(n, [](int i, int j, double d) { return i == j ? 1000.0 + d : d; });
    else if (target == "bad")
        generate(n, [](int i, int j, double d) {return d;});
    else if (target == "Gilbert")
        generate(n, [](int i, int j, double d) {return static_cast<double>(1. / (1 + i + j));});
    else cout << "undefined matrix" << endl;
}
