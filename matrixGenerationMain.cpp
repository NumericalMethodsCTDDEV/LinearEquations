#include <bits/stdc++.h>

using namespace std;

const int N = 100;

int main(int argc, char *argv[])
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1, N);

    double dom = 0;
    if (argc == 1) dom = 0;
    else
    {
        dom = 1000;
    }

    size_t n = 10;

    cout << n << endl;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            double ok = 0.0;
            if (i == j) ok = dom;
            cout << dis(gen) + ok << " ";
        }
        cout << endl;
    }
    for (size_t i = 0; i < n; ++i)
        cout << dis(gen) << endl;
}
