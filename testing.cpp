#include "gtest/gtest.h"
#include "linearSystemsSolver.h"
#include <limits.h>

using namespace linearSystemsSolver;

TEST(testing, simple_test1)
{
    matrix_t a = { {0, 0, 0}, {0, 0, 0}};
    answer_t response = solve(a, "gauss");
    EXPECT_TRUE(response.amountOfSolutions == linearSystemsSolver::Inf);
}

TEST(testing, simple_test2)
{
    matrix_t a = { {0, 0, 1}, {0, 0, 1}};
    answer_t response = solve(a, "gauss");
    EXPECT_TRUE(response.amountOfSolutions == 0);
}

TEST(testing, simple_test3)
{
    matrix_t a = { {3, 1, 2}, {1, -2, 3}};
    answer_t response = solve(a, "jacobi");
    EXPECT_TRUE(response.amountOfSolutions == 1);
    std::vector<double> right_ans = {1, -1};
    EXPECT_TRUE(response.solution == right_ans);
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
