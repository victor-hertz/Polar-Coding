#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;
template<typename T>

vector<size_t> tag_sort(const vector<T>& v)
{
    vector<size_t> result(v.size());
    iota(begin(result), end(result), 0);
    sort(begin(result), end(result),
            [&v](const auto & lhs, const auto & rhs)
            {
                return v[lhs] < v[rhs];
            }
    );
    return result;
}

const vector<double> ZERO_MESSAGE = {1, 0, 0};
const vector<double> ONE_MESSAGE = {0, 1, 0};
const vector<double> UN_MESSAGE = {0.5, 0.5, 0};

bool random_bool();
vector<bool> random_bool_list(int size);
vector<bool> zero_bool_list(int size);
vector<int> random_integer_list(int max, int size);
int reverse_bits(int n, int m);
vector<bool> reverse_bits(vector<bool> &l, int m);
int nb_differences(const vector<bool> &a, const vector<bool> &b);
vector<int> pos_difference(const vector<bool> &a, const vector<bool> &b);

vector<vector<double>> print_stats(vector<vector<double>> lik);
void print_comp(vector<vector<double>> lik_std, vector<vector<double>> lik_par, vector<int> pos_error);

void print_list(vector<bool> l);
void print_list(vector<double> l);
void print_list(vector<int> l);
void print_list(vector<unsigned long> l);
void print_list(vector<vector<double>> l);
void print_list(vector<vector<int>> l);

#endif // UTILS_H
