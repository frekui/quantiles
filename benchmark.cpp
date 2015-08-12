// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

// Code used for benchmarking SbQuantiles.

#define SBQ_TEST 0

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

#include <math.h>

#include "ExactQuantiles.h"
#include "SbQuantiles.h"
#include "TreeNode.h"

using namespace std;

class Timer
{
public:
    Timer() : start(std::chrono::system_clock::now())
    { }

    double get_elapsed() const
    {
        auto end = chrono::system_clock::now();
        chrono::duration<double> dur(end - start);
        return dur.count();
    }
    
private:
    std::chrono::time_point<std::chrono::system_clock> start;
};

enum class Compress {none, all, automatic};

inline std::ostream& operator<< (std::ostream& o, const Compress& c)
{
    switch (c) {
        case Compress::none: return o << "none";
        case Compress::all: return o << "all";
        case Compress::automatic: return o << "automatic";
        default: assert(false);
    }
}

void benchmark(int UniverseBits, Compress compress, int seed = 0)
{
    if (SBQ_TEST)
        abort();
    
    cout << "Benchmark,"
        << " UniverseBits: " << UniverseBits
        << " compress: " << compress
        << " seed: " << seed
        << endl;
    SbQuantiles h(UniverseBits, 0.01);
    h.print_params();
    cout << '\n';
    uint32_t U_size = TreeNode::u_size(UniverseBits);
    h.set_auto_compress(compress == Compress::automatic ? true : false);
    std::mt19937 gen;
    gen.seed(seed);
    std::uniform_int_distribution<> dis(0, U_size-1);
    int items = 500000;
    
    Timer timer;
    for (int i = 0; i < items; i++) {
        int x = dis(gen);
        h.insert(x);
        
        if (compress == Compress::all)
            h.compress_from_root();
    }
    
    printf("Inserted %d items, %.1fs. %.1f items/s compress: %s\n",
           items,
           timer.get_elapsed(),
           items/timer.get_elapsed(),
           to_string(compress).c_str());
    Timer timer_compress;
    h.compress_from_root();
    printf("Final compress took %.1f s\n", timer_compress.get_elapsed());
    cout << "Size: " << h.size() << " N: " << h.get_N() << endl;
    cout << endl;
}

int main(int, const char**)
{
    benchmark(30, Compress::automatic, 1);
    return 0;
}
