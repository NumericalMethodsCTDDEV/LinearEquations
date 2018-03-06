#include <vector>
#include <cassert>
#include "linearSystemsSolver.h"

#include <iostream>
#include <cmath>

using namespace std;

typedef std::vector<double> v_t;
namespace matrix_utils {
    double product(const v_t &a, const v_t &b) {
        assert(a.size() == b.size());
        double res = 0;
        for (size_t w = 0; w < a.size(); w++) {
            res += a[w] * b[w];
        }
        return res;
    }
    v_t product(const vector<v_t> &a, const v_t &b) {
        assert(a.size() == b.size());
        v_t res(a.size(), 0);
        for (size_t e = 0; e < a.size(); e++) {
            for (size_t r = 0; r < a.size(); r++) {
                res[e] += b[r] * a[e][r];
            }
        }
        return res;
    }
    vector<v_t> product(const vector<v_t> &a, const vector<v_t> &b) {
        assert(a.size() == b.size());
        size_t size = a.size();
        vector<v_t> res(size, v_t(size, 0));
        
        for (size_t w = 0; w < size; w++) {
            for (size_t e = 0; e < size; e++) {
                for (size_t r = 0; r < size; r++) {
                    res[w][r] += a[w][e] * b[e][r];
                }
            }
        }
        return res;
    }
    v_t mult_and_add(double coeff_1, const v_t &a, double coeff_2, const v_t &b) {
        assert(a.size() == b.size());
        v_t res(a.size(), 0);
        for (size_t w = 0; w < a.size(); w++) {
            res[w] = a[w] * coeff_1 + b[w] * coeff_2;
        }
        return res;
    }
    vector<v_t> transpose(const vector<v_t> &a) {
        vector<v_t> res(a.size(), v_t(a.size(), 0));
        
        for (size_t e = 0; e < a.size(); e++) {
            for (size_t r = 0; r < a.size(); r++) {
                res[e][r] = a[r][e];
            }
        }
        return res;
    }
}
using namespace matrix_utils;

namespace {
    const double EPS = 1e-9;
    
    void symmetrization(vector<v_t> &a, v_t &b) {
        vector<v_t> a_t = transpose(a);
        a = product(a_t, a);
        b = product(a_t, b);
    }
    v_t calc_delta(const vector<v_t> &a, const v_t &b, const v_t &x) {
        return mult_and_add(1, product(a, x), -1, b);
    }
    double calc_beta(const vector<v_t> &a, const v_t &dir, const v_t &new_delta) {
        v_t mult = product(a, dir);
        return product(mult, new_delta) / product(mult, dir);
    }
    double calc_alpha(const vector<v_t> &a, const v_t &dir, const v_t &delta, const v_t &x, bool &error) {
        double A = product(product(a, dir), dir) / 2;
        double B = product(delta, dir);
        
        if (A <= 0) {
//            cout << A << " !!!\n";
            error = true;
            return 0;
        }
        return -B / (2 * A);
    }
    double vector_norm(const v_t &a) {
        double res = 0;
        for (double val : a) {
            res += abs(val);
        }
        return res;
    }
    double matrix_norm(const vector<v_t> &a) {
        double res = -1;
        
        for (size_t e = 0; e < a.size(); e++) {
            double sum = 0;
            for (size_t r = 0; r < a.size(); r++) {
                sum += abs(a[r][e]);
            }
            res = max(res, sum);
        }
        return res;
    }
    double real_error(const v_t &solution, const vector<v_t> &original_a, const v_t &original_b) {
        v_t delta = mult_and_add(1, product(original_a, solution), -1, original_b);
        
        cout << vector_norm(delta) << "\n";
        cout << vector_norm(delta) / vector_norm(solution) << "\n";
        
        return min(vector_norm(delta), vector_norm(delta) / vector_norm(solution));
    }
    int check_solution(const v_t &solution, const vector<v_t> &original_a, const v_t &original_b) {
        if (real_error(solution, original_a, original_b) <= EPS) {
            return 1;
        } else {
            return 3;
        }
    }
}

int descent(vector<v_t> a, v_t &ans) {
    assert(a.size() == a.front().size() - 1);
    v_t b;
    for (v_t &vect : a) {
        b.push_back(vect.back());
        vect.pop_back();
    }
    vector<v_t> original_a = a;
    v_t original_b = b;
    
//    double corrected_eps = EPS / matrix_norm(transpose(a));
//    cout << corrected_eps << "\n";
    
    symmetrization(a, b);
    
    v_t x(a.size(), 0);
    for (int iteration = 0; iteration < 5; iteration++) {
        v_t delta = calc_delta(a, b, x);
        v_t dir = delta;
        
        for (double &w : dir) {
            w = -w;
        }
        if (vector_norm(delta) < EPS) {
            ans = x;
            return 1;
        }
        
        for (size_t w = 0; w < a.size(); w++) {
            bool error = false;
            double alpha = calc_alpha(a, dir, delta, x, error);
            
            if (error) {
                return 3;
            }
            v_t new_x = mult_and_add(1, x, alpha, dir);
            v_t new_delta = calc_delta(a, b, new_x);
            
            double n_2 = vector_norm(new_x);
            double n_3 = vector_norm(new_delta);
            
            if (((n_3 < EPS) || (n_3 / n_2 < EPS))) {
                ans = new_x;
                return 1;
            }
            
            v_t new_dir = mult_and_add(-1, new_delta, calc_beta(a, dir, new_delta), dir);
            
            x     = new_x;
            delta = new_delta;
            dir   = new_dir;
        }
    }
    return 2;
}
