// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

// Tests for SbQuantiles and related classes.

#if !SBQ_TEST
#  error "SBG_DEBUG must be 1"
#endif

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

#include "ExactQuantiles.h"
#include "SbQuantiles.h"
#include "TreeNode.h"

using namespace std;

static const bool verbose = false;

enum class Compress {none, all, automatic};

static inline std::ostream& operator<< (std::ostream& o, const Compress& c)
{
    switch (c) {
        case Compress::none: return o << "none";
        case Compress::all: return o << "all";
        case Compress::automatic: return o << "automatic";
        default: assert(false);
    }
}

static void random_test(int UniverseBits, Compress compress, int seed = 0)
{
    uint32_t U_size = TreeNode::u_size(UniverseBits);
    constexpr double epsilon = 0.5;
    SbQuantiles h(UniverseBits, epsilon);
    cout << "Running test with random data. compress: " << compress << ' ';
    h.print_params();
    cout << "... " << flush;
    if (compress != Compress::automatic)
        h.set_auto_compress(false);
    ExactQuantiles e;
    std::mt19937 gen;
    gen.seed(seed);
    std::uniform_int_distribution<> dis(0, U_size-1);
    // bool print_tree = U_size <= 32;
    bool print_tree = false;
    
    for (int i = 0; i < 1000; i++) {
        int x = dis(gen);
        if (verbose)
            cout << i << "\tinserting " << x << endl;
        e.insert(x);
        h.insert(x);
        assert(h.test_invariants());
        
        if (compress == Compress::all) {
            if (verbose) {
                cout << "Compressing..." << endl;
                h.dump();
            }
            if (print_tree)
                h.print_tree();
            h.compress_from_root();
            assert(h.test_invariants());
            assert(h.test_accuracy(e, seed, 1000, 100));
            if (verbose)
                cout << "Size: " << h.size() << " N: " << h.get_N() << endl;
        }
    }
    
    if (print_tree)
        h.print_tree();
    cout << "ok" << endl;
}

static void do_tests()
{
    bool simple_tests = true;
    int test_iterations = simple_tests ? 10 : 1;
    
    for (int i = 0; i < test_iterations; i++)
        random_test(3, Compress::all, i+3);
    random_test(4, Compress::none);
    random_test(4, Compress::all);
    random_test(5, Compress::none);
    for (int i = 0; i < test_iterations; i++)
        random_test(4, Compress::all, i+3);
    random_test(5, Compress::all);
    for (int i = 0; i < test_iterations; i++)
        random_test(5, Compress::all, i+3);

    // Example from the paper.
    int UniverseBits = 4;
    uint32_t U_size = TreeNode::u_size(UniverseBits);
    SbQuantiles h(UniverseBits, 0.5);
    cout << "Running test with example from paper ";
    h.print_params();
    cout << "... " << flush;
    h.set_auto_compress(false);
    ExactQuantiles e;
    for (int x : {1,2,2,2,3,4,5,6,6,6,10,12,14,14,15,16}) {
        // We use [0,15] as the universe instead of [1,16]
        x = x - 1;
        e.insert(x);
        if (verbose)
            cout << "Inserting " << x << endl;
        h.insert(x);
        if (verbose)
            h.dump();
        assert(h.test_invariants());
    }
    
    if (verbose)
        h.print_tree();
    h.compress_from_root();
    if (verbose)
        h.print_tree();
    assert(h.test_invariants());
    cout << "ok." << endl;
    if (verbose) {
        for (int x = 0; x < U_size; x++) {
            cout << "x: " << x
                 << "\texact rank: " << e.rank(x)
                 << "\tapprox rank: " << h.rank(x)
                 << "\tdiff: " << e.rank(x) - h.rank(x)
                 << endl;
        }
    }
}

// A test case which shows that the pseudo code in CompressTree in Figure 3
// in the paper has a bug. When the code from the paper is used the assert
// statement at the end of the function fails because the L value for node 12
// in the tree is computed as 6, but it should be 8.
//
// The problem has been acknowledged by Graham Cormode, one of the co-authors
// of the paper.
static void bug_small()
{
    SbQuantiles h(3, 0.5);
    cout << "Running test bug_small ";
    h.print_params();
    cout << "... " << flush;
    h.set_auto_compress(false);
    h.insert(0);
    h.insert(0);
    h.insert(0);
    h.insert(0);
    h.insert(0);
    h.insert(0);
    h.insert(5);
    h.insert(5);
    if (verbose) {
        cout << "After 6x0, 2x5\n";
        h.print_tree();
    }
    h.compress_from_root();
    if (verbose) {
        cout << "After 6x0, 2x5, compress\n";
        h.print_tree();
    }
    assert(h.test_invariants());
    cout << "ok." << endl;
}

int main(int, const char**)
{
    test_assert_enable();
    bug_small();
    testTreeNode();
    do_tests();
    return 0;
}
