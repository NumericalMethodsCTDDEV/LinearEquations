#include "gtest/gtest.h"
#include "linearSystemsSolver.h"
#include <limits.h>
#include <cmath>

using namespace linearSystemsSolver;

TEST(testing, test_inf)
{
    std::vector<std::string> allMethods = {{"gauss"}};
    matrix_t a = { {0, 0, 0}, {0, 0, 0}};
    std::vector<answer_t> responses;
    for (const auto &name : allMethods)
        responses.push_back(solve(a, name.c_str()));
    for (const auto &ri : responses)
        EXPECT_TRUE(ri.amountOfSolutions == Inf);
}

TEST(testing, test_zero)
{
    std::vector<std::string> allMethods = {{"gauss"}};
    matrix_t a = { {0, 0, 1}, {0, 0, 1}};
    std::vector<answer_t> responses;
    for (const auto &name : allMethods)
        responses.push_back(solve(a, name.c_str()));
    for (const auto &ri : responses)
        EXPECT_TRUE(ri.amountOfSolutions == 0);
}

TEST(testing, just_test)
{
    std::vector<std::string> allMethods = {{"gauss"}};
    matrix_t a = { {1, 1, 2}, {1, -1, 0}};
    std::vector<double> right_ans = {1, 1};
    std::vector<answer_t> responses;
    for (const auto &name : allMethods)
        responses.push_back(solve(a, name.c_str()));
    for (const auto &ri : responses)
    {
        EXPECT_TRUE(ri.amountOfSolutions == 1);
        EXPECT_TRUE(ri.solution == right_ans);
    }
}


TEST(testing, just_test2)
{
    std::vector<std::string> allMethods = {{"gauss"}, {"seidel"}, {"sequantialRelaxation"}};
    matrix_t a = { {3, 1, 2}, {1, -2, 3}};
    std::vector<double> right_ans = {1, -1};
    std::vector<answer_t> responses;
    for (const auto &name : allMethods)
        responses.push_back(solve(a, name.c_str()));
    for (const auto &ri : responses)
    {
        EXPECT_TRUE(ri.amountOfSolutions == 1);
        for (size_t i = 0; i < ri.solution.size(); ++i)
            if (ri.solution[i] < right_ans[i] - EPS || ri.solution[i] > right_ans[i] + EPS)
                EXPECT_TRUE(false);
    }
}

TEST(testing, jacobi_test)
{
    matrix_t a = { {3, 1, 2}, {1, -2, 3}};
    answer_t response = solve(a, "jacobi");
    EXPECT_TRUE(response.amountOfSolutions == 1);
    std::vector<double> right_ans = {1, -1};
    if (response.solution[0] < right_ans[0] - EPS || response.solution[0] > right_ans[0] + EPS)
        EXPECT_TRUE(false);
}


TEST(testing, seidel_test)
{
    matrix_t a = { {3, 1, 2}, {1, -2, 3}};
    answer_t response = solve(a, "seidel");
    EXPECT_TRUE(response.amountOfSolutions == 1);
    std::vector<double> right_ans = {1, -1};
    if (response.solution[0] < right_ans[0] - EPS || response.solution[0] > right_ans[0] + EPS)
        EXPECT_TRUE(false);
}

TEST(testing, seq_relaxation_test)
{
    matrix_t a = { {3, 1, 2}, {1, -2, 3}};
    answer_t response = solve(a, "sequantialRelaxation");
    EXPECT_TRUE(response.amountOfSolutions == 1);
    std::vector<double> right_ans = {1, -1};
    if (response.solution[0] < right_ans[0] - EPS || response.solution[0] > right_ans[0] + EPS)
        EXPECT_TRUE(false);
}

namespace
{
    const size_t RANDOM_TESTS_SIZE = 10000;
    void printSystem(const matrix_t &a)
    {
        for (auto &vi : a)
        {
            for (auto i : vi)
                std::cerr << i << " ";
            std::cerr << std::endl;
        }
        std::cerr << std::endl;
    }
}

TEST(testing, multi_random_test)
{
    std::vector<std::string> allMethods = getAllAvailableMethods();
    srand(time(0));
    size_t cntGood = 0;
    for (size_t t = 0; t < RANDOM_TESTS_SIZE; ++t)
    {
        size_t n = rand() % 10 + 1;
        std::vector<std::vector<double>> a(n);
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                a[i].push_back(double(rand() % 10 + (i == j ? 100 : 1)));
            }
            a[i].push_back(double(rand() % 10 + 1));
        }
        //        printSystem(a);
        std::vector<answer_t> responses;
        for (const auto &name : allMethods)
            responses.push_back(solve(a, name.c_str()));
        for (const auto &ri : responses)
        {
            for (const auto &rj : responses)
            {
                //                EXPECT_TRUE(ri.amountOfSolutions == rj.amountOfSolutions);
                if (ri.amountOfSolutions == 1 && rj.amountOfSolutions == 1)
                {
                    ++cntGood;
                    for (size_t k = 0; k < ri.solution.size(); ++k)
                        EXPECT_TRUE(std::fabs(ri.solution[k] - rj.solution[k]) < 2 * EPS);
                }
            }
        }
    }
    //    std::cerr << cntGood / (allMethods.size() * allMethods.size()) << std::endl;
}

TEST(testing, transpose_test)
{
    matrix_t a = { {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
    matrix_t ta = { {1, 1, 1, 1, 1}, {2, 2, 2, 2, 2}, {3, 3, 3, 3, 3}, {4, 4, 4, 4, 4}, {5, 5, 5, 5, 5}};
    matrix_t my_ta = transpose(a);
    EXPECT_TRUE(my_ta == ta);
    matrix_t my_a = transpose(ta);
    EXPECT_TRUE(my_a == a);
    a = { {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
    ta = { {1, 1, 1, 1}, {2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};
    my_ta = transpose(a);
    EXPECT_TRUE(ta == my_ta);
    a = {{1}};
    my_ta = transpose(a);
    EXPECT_TRUE(a == my_ta);
}

TEST(testing, determinant_test)
{
    matrix_t a = {{1, 2}, {2, 1}};
    EXPECT_TRUE(determinant(a) == -3);
    a = {{0, 1}, {0, 1}};
    EXPECT_TRUE(determinant(a) == 0);
    a = {{10.0001}};
    EXPECT_TRUE(determinant(a) == 10.0001);
    a = { {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
    EXPECT_TRUE(determinant(a) == 0);
    a = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 4}, {1, 2, 3, 3, 3}, {1, 2, 2, 2, 2}, {5, 4, 3, 2, 1}};
    EXPECT_TRUE(determinant(a) == 6);
}

TEST(testing, descent_test) {
    matrix_t a = {
        {6.000000, 7.000000, 3.000000, -5.000000, -2.000000, 5.000000, 5.000000, 9.000000, 4.000000, -8.000000, 2.000000},
        {-5.000000, -6.000000, -3.000000, 5.000000, 5.000000, -2.000000, -3.000000, 5.000000, 8.000000, 9.000000, 5.000000},
        {-1.000000, 2.000000, -6.000000, 2.000000, 8.000000, 2.000000, -8.000000, -4.000000, -10.000000, 8.000000, 0.000000},
        {-10.000000, -1.000000, 8.000000, -8.000000, 4.000000, -9.000000, 2.000000, -1.000000, -1.000000, 2.000000, 7.000000},
        {-10.000000, 0.000000, -8.000000, 7.000000, -10.000000, -8.000000, -3.000000, 8.000000, 3.000000, -4.000000, 10.000000},
        {-8.000000, -2.000000, -10.000000, -8.000000, 1.000000, -7.000000, 6.000000, -7.000000, 9.000000, 2.000000, -3.000000},
        {1.000000, 9.000000, 0.000000, 4.000000, 0.000000, -3.000000, -10.000000, 6.000000, -5.000000, -10.000000, 8.000000},
        {-2.000000, 2.000000, -8.000000, -5.000000, 10.000000, -4.000000, 0.000000, 7.000000, -2.000000, 9.000000, -3.000000},
        {-2.000000, 7.000000, 6.000000, 8.000000, 4.000000, 4.000000, -7.000000, 10.000000, -6.000000, 0.000000, 6.000000},
        {4.000000, 7.000000, 6.000000, -10.000000, -6.000000, -4.000000, -5.000000, 2.000000, 7.000000, 6.000000, -1.000000}
    };
    std::vector<double> res;
    EXPECT_TRUE(determinant(a) != 0);
    EXPECT_TRUE(solve(a, "descent").amountOfSolutions == 1);
}
