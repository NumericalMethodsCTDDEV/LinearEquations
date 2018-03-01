#include <gtest/gtest.h>
#include "gauss.h"
#include <limits.h>

TEST(testing, simple_test1)
{
    std::vector<std::vector<double>> a = { {0, 0, 0}, {0, 0, 0}};
    std::vector<double> ans;
    int status = gauss(a, ans);
    EXPECT_TRUE(status == INT_MAX - 100);
}

TEST(testing, simple_test2)
{
    std::vector<std::vector<double>> a = { {0, 0, 1}, {0, 0, 1}};
    std::vector<double> ans;
    int status = gauss(a, ans);
    EXPECT_TRUE(status == 0);
}

TEST(testing, simple_test3)
{
    std::vector<std::vector<double>> a = { {1, 1, 2}, {1, -1, 0}};
    std::vector<double> ans;
    int status = gauss(a, ans);
    EXPECT_TRUE(status == 1);
    std::vector<double> right_ans = {1, 1};
    EXPECT_TRUE(right_ans == ans);
}
