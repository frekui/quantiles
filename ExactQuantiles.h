// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#ifndef ExactQuantiles_h
#define ExactQuantiles_h

#include "common.h"
#include <map>

// A small class to exactly compute the rank of inserted items. It is
// slow and a memory hog, but the implementation is small and simple
// and it is useful when testing SbQuantiles.
class ExactQuantiles
{
    using U = uint32_t;
    std::map<U, unsigned> counts;

public:
    void insert(U x);
    int rank(uint64_t x) const;
    int get_N() const;
    U quantile(double p) const;
};

#endif
