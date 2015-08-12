// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

// A small test program which demonstrates how SbQuantiles can be
// used.

#include "SbQuantiles.h"
#include <random>
#include <cstdio>

using namespace std;

void print_quantile(SbQuantiles h, int U_size, double phi)
{
    // Get the approximate phi quantile.
    int x_hat = h.quantile(phi);
    
    // As we are using an uniform distribution we can compute the correct
    // answer (modulo some random noise as we are picking values at random).
    int x = phi*U_size;
    
    int rank_x_hat = h.get_N()*(double(x_hat)/U_size);
    int rank_x = h.get_N()*(double(x)/U_size); // can also be written as h.get_N()*phi.

    // Note: No bounds are given on |x_hat - x|. The accuracy is in terms of
    // rank and is given by the inequalities
    // true_rank(x)*(1 − epsilon) <= phi*N <= true_rank(x + 1)*(1 + epsilon).
    printf("Approx %.2f quantile x_hat = %10d. rank(x) = %6d, "
           "rank(x_hat) = %6d diff in rank: %d\n",
           phi, x_hat, rank_x, rank_x_hat, rank_x - rank_x_hat);
}

int main(void)
{
    // Size of our universe. The items we want to insert must lie in the range
    // [0, 2**UniverseBits).
    int UniverseBits = 30;
    int U_size = TreeNode::u_size(UniverseBits);
    
    // Initialize the data structure. The second argument is the accuracy
    // parameter epsilon. Lower epsilon means higher accuracy but also a larger
    // and slower data structure.
    double epsilon = 0.01;
    SbQuantiles h(UniverseBits, epsilon);
    cout << "Using U = 2**" << UniverseBits
         << " and epsilon = " << epsilon
         << endl;

    mt19937 gen;
    random_device rd;
    gen.seed(rd());
    uniform_int_distribution<> dis(0, U_size-1);
    
    const int items = 500000;
    cout << "Will insert " << items << " uniformly distributed integers.\n";
    for (int i = 0; i < items; i++)
        h.insert(dis(gen));

    // Since the data structure is designed for fully biased rank queries the
    // accuracy will not be the same for different quantiles. In this
    // implementation the accuracy will go up when lower quantiles are asked
    // for. In particular, the accuracy of the last two quantiles below will
    // not be great compared to the others.
    for (double p : {0.0, 0.01, 0.05, 0.5, 0.95, 0.99})
        print_quantile(h, U_size, p);

    return 0;
}
