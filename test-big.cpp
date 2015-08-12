// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

// Test program which tests with a large universe. This program does not call
// test_assert_enable as it would become too slow if TEST_ASSERT were enabled.

#include "common.h"

#include <iostream>
#include <map>
#include <set>
#include <numeric>
#include <random>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <cstdio>

#include <assert.h>
#include <math.h>

#include "SbQuantiles.h"
#include "ExactQuantiles.h"

using namespace std;

// Print histogram on stdout.
// FIXME: Move to example.cpp?
template<class Hist>
static void print_histogram(const Hist& h, int buckets)
{
    uint64_t start = h.quantile(0.01);
    uint64_t end = h.quantile(0.99);
    uint64_t delta = (end-start)/buckets;
    if (delta == 0)
        delta = 1;
    if (start < delta)
        start = delta;

    const int num_stars = 100;
    int prev_rank = 0;
    printf("%-7s %-10s\n", "Bucket", "#items");
    for (int i = 0; start < end; start += delta, i++) {
        int cur_rank = h.rank(start);
        string stars(num_stars*double(cur_rank - prev_rank)/h.get_N(), '*');
        printf("%-7d %-10d %-4.1f%% %s\n",
               i,
               cur_rank - prev_rank,
               100*double(cur_rank-prev_rank)/h.get_N(),
               stars.c_str());
        prev_rank = cur_rank;
    }
}


template<class Dist>
static void random_test(int UniverseBits, double epsilon,
                        int seed, Dist dist)
{
    const int items_to_insert = 1000000;
    const int items_to_test = 100;
    const int quantiles_to_test = 10;

    uint32_t U_size = TreeNode::u_size(UniverseBits);
    SbQuantiles h(UniverseBits, epsilon);
    h.print_params();
    cout << '\n';
    ExactQuantiles e;

    mt19937 gen;
    gen.seed(seed);

    cout << "Inserting " << items_to_insert << " items." << endl;
    for (int i = 0; i < items_to_insert; i++) {
        int x = dist(gen);
        if (x >= U_size)
            continue;
        e.insert(x);
        h.insert(x);
    }
    cout << "Number of tuples: " << h.size() << " (after insertions)." << endl;

    cout << "Testing accuracy guarantees... " << flush;
    assert(h.test_accuracy(e, seed, items_to_test, quantiles_to_test));
    cout << "ok." << endl;

    cout << "\nApproximate histogram\n";
    print_histogram(h, 20);
    cout << "\nExact histogram\n";
    print_histogram(e, 20);
}

int main(void)
{
    mt19937 gen;
    gen.seed(123);
    int UniverseBits = 31;
    uint32_t U_size = TreeNode::u_size(UniverseBits);

    cout << "Testing with uniform input distribution" << endl;
    random_test(UniverseBits, 0.01, 123,
                uniform_int_distribution<uint32_t>(0, U_size-1));

    cout << "\nTesting with geometric input distribution" << endl;
    random_test(UniverseBits, 0.01, 654,
                geometric_distribution<uint32_t>(0.005));

    return 0;
}
