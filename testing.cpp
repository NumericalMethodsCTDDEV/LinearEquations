#include <gtest/gtest.h>
#include "linearSystemsSolver.h"
#include <limits.h>

TEST(testing, simple_test1)
{
    linearSystemsSolver solver;
    system_t a = { {0, 0, 0}, {0, 0, 0}};
    answer_t response = solver.solve(a, "gauss");
    EXPECT_TRUE(response.amountOfSolutions == INT_MAX - 100);
}

TEST(testing, simple_test2)
{
    linearSystemsSolver solver;
    system_t a = { {0, 0, 1}, {0, 0, 1}};
    answer_t response = solver.solve(a, "gauss");
    EXPECT_TRUE(response.amountOfSolutions == 0);
}

TEST(testing, simple_test3)
{
    linearSystemsSolver solver;
    system_t a = { {1, 1, 2}, {1, -1, 0}};
    answer_t response = solver.solve(a, "gauss");
    EXPECT_TRUE(response.amountOfSolutions == 1);
    std::vector<double> right_ans = {1, 1};
    EXPECT_TRUE(response.solution == right_ans);
}
