// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#include <assert.h>

#include "ExactQuantiles.h"

using namespace std;

void ExactQuantiles::insert(U x)
{
    auto it = counts.find(x);
    if (it == counts.end())
        counts[x] = 1;
    else
        it->second += 1;
}

int ExactQuantiles::rank(uint64_t x) const
{
    int sum = 0;
    for (auto p : counts) {
        if (p.first < x)
            sum += p.second;
        else
            break;
    }
    
    return sum;
}

int ExactQuantiles::get_N() const
{
    int N = 0;
    for (auto p : counts)
        N += p.second;

    return N;
}

ExactQuantiles::U ExactQuantiles::quantile(double phi) const
{
    assert(0 <= phi && phi <= 1);
    uint32_t prev = 0;
    unsigned sum = 0;
    unsigned N = get_N();
    for (auto p : counts) {
        sum += p.second;
        if (sum > phi*N)
            break;
        prev = p.first;
    }

    return prev;
}
