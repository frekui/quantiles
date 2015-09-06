# What is this?
This repository contains a C++11 implementation of the data structure
described in Section 4.1 of the paper

1. G. Cormode, F. Korn, S. Muthukrishnan, and D. Srivastava. Space- and time-efficient deterministic algorithms for biased quantiles over data streams. In ACM Principles of Database Systems (PODS), 2006.

The intended use-case is that we have a large stream of integers and we want to efficiently approximate rank and quantile queries for the stream without storing the whole stream.

The data structure supports the following operations:

## construction(b, &epsilon;)
Preconditions: 0 < b &le; 31 and 0 < &epsilon; &le; 0.5.

b is the number of bits in the universe and U = 2<sup>b</sup> is the size of the universe. In this implementation U must be less
than 2<sup>31</sup> (i.e., 31 bit integers).

&epsilon; is a parameter which controls the trade-off between space
usage and accuracy. The space used by the data structure is bounded by
O(log(U) &middot; 1/&epsilon; &middot; (log(&epsilon; N) + log(U))). See rank and quantile below for the accuracy guarantees.

## insert(x)
Precondition: x is an integer such that 0 &le; x < U.

Insert x. Amortized time complexity: O(log(U)).

## rank(x)
Precondition: x is an integer such that 0 &le; x < U.

Compute the approximate rank of x. Let
true\_rank(x) = number of inserted items y such that y < x.
Then it is guaranteed that abs(rank(x) - true\_rank(x)) < &epsilon; &middot; true\_rank(x).

## quantile(p)
Precondition: 0 &le; p &le; 1.

Compute an approximate p quantile. It is guaranteed that 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;true\_rank(x) &middot; (1 − &epsilon;) &le; pN &le; true\_rank(x + 1) &middot; (1 + &epsilon;)

where x = quantile(p) and N is the number of items that have been inserted. quantile(p) is implemented by binary searching for a suitable x using rank(&middot;).

# Structure of the code
The code for the data structure is contained in [`SbQuantiles.h`](SbQuantiles.h),
[`SbQuantiles.cpp`](SbQuantiles.cpp), and [`TreeNode.h`](TreeNode.h).

Any references to theorems, lemmas, definitions, equations, figures, and
sections in comments in the code refers to the corresponding piece of
text in [1].

## [`example.cpp`](example.cpp)
A fairly small example which show how SbQuantiles can be
used.

## [`benchmark.cpp`](benchmark.cpp)
Used for benchmarking.

## [`test.cpp`](test.cpp) and [`test-big.cpp`](test-big.cpp)
Tests.

# License

The MIT license, see [`LICENSE.txt`](LICENSE.txt) for details.
Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.
