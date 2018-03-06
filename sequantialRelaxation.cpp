#include <vector>
#include <iostream>
#include <cmath>
#include "linearSystemsSolver.h"

const double W = 1; // коэффициент релаксации:
                    // 0 < W < 1 - последовательная нижняя релаксация
                    // W = 1 - обычный метод Зайделя
                    // 1 < W < 2 - последовательная верхняя релаксация

bool sequantialRelaxationSatisfyEst(double est, const std::vector<double>& x, const std::vector<double>& xNext) {
    double norm = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        norm += (xNext[i] - x[i]) * (xNext[i] - x[i]);
    }
    norm = sqrt(norm);
    return norm <= est;
}

int sequantialRelaxation(std::vector<std::vector<double>> a, std::vector<double> &ans) {

    size_t n = a.size();

    // находим матрицы nxn b, b1, b2 и свободный член c
    std::vector<std::vector<double>> b(n);
    std::vector<std::vector<double>> b1(n); // нижнетреугольная
    std::vector<std::vector<double>> b2(n); // верхнетреугольная
    std::vector<double> c;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            b[i].push_back(i != j ? -a[i][j] / a[i][i] : 0);
            b1[i].push_back(i > j ? b[i][j] : 0);
            b2[i].push_back(i < j ? b[i][j] : 0);
        }
        c.push_back(a[i].back() / a[i][i]);
    }

    // проверяем достаточное условие сходимости
    bool convergance = linearSystemsSolver::norm(b1) + linearSystemsSolver::norm(b2) < 1;

    // находим апостериорную оценку сходимости
    double est = std::min(linearSystemsSolver::EPS, (1 - linearSystemsSolver::norm(b)) / linearSystemsSolver::norm(b2) * linearSystemsSolver::EPS);

    // задаем начальный вектор
    std::vector<double> x, xNext;
    for (size_t i = 0; i < n; ++i) {
        xNext.push_back(0);
    }

    // находим следующее приближение, пока оно и текущее не будут удовлетворять апостериорной оценке
    bool begin = true;
    std::vector<double> xNextEst;
    size_t iterations = linearSystemsSolver::ITERATIONS;
    while (begin || (convergance ? !sequantialRelaxationSatisfyEst(est, x, xNext) : iterations > 0)) {
        begin = false;
        --iterations;
        x = xNext;
        xNext.clear();
        xNextEst.clear();
        for (size_t i = 0; i < n; ++i) {
            double xNextEstElement = c[i];
            for (size_t j = 0; j < n; ++j) {
                xNextEstElement += b[i][j] * (i > j ? xNext[j] : x[j]);
            }
            xNextEst.push_back(xNextEstElement);
        }
        for (size_t i = 0; i < n; ++i) {
            xNext.push_back(W * xNextEst[i] + (1 - W) * x[i]);
        }
    }

    // получаем ответ
    ans = xNext;
    if (convergance) {
        return 1;
    }
    return 2;
}