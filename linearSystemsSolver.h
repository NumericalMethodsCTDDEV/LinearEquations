#ifndef LINEARSYSTEMSSOLVER_H
#define LINEARSYSTEMSSOLVER_H

#include <string>
#include <vector>
#include <limits.h>
#include <unordered_map>

struct answer_t
{
    std::string status;
    std::vector<double> solution;
    int amountOfSolutions;
};

using system_t = std::vector<std::vector<double>>;

class linearSystemsSolver
{
    static int gauss(std::vector < std::vector<double> > a, std::vector<double> &ans);
    static const int Inf = INT_MAX - 100;

    std::string parseStatus(int st)
    {
        switch (st)
        {
        case 0:
            return "No solutions";
        case 1:
            return "One solution";
        case Inf:
            return "Inf solutions";
        default:
            return "Something went wrong";
        }
    }

    using solveMethdo_t = int (system_t, std::vector<double> &ans);
    std::unordered_map<std::string, solveMethdo_t *> nameToMethod;

public:

    linearSystemsSolver()
    {
        nameToMethod["gauss"] = &linearSystemsSolver::gauss;
    }

    answer_t solve(system_t &system, const char *methodName)
    {
        std::string name(methodName);
        std::vector<double> ans;
        if (!nameToMethod.count(name))
            return {parseStatus(-1), {}, 0};
        solveMethdo_t *f = nameToMethod[name];
        int amount = f(system, ans);
        std::string status = parseStatus(amount);
        return {status, ans, amount};
    }

};

#endif // LINEARSYSTEMSSOLVER_H
