// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#ifndef SbQuantiles_h
#define SbQuantiles_h

#include "common.h"

#include <iostream>
#include <unordered_map>
#include <numeric>
#include <random>
#include <algorithm>

#include <assert.h>
#include <math.h>

#include "TreeNode.h"

class ExactQuantiles;
// This class is an implementation of the data structure described in Section
// 4.1 of the paper
//
// G. Cormode, F. Korn, S. Muthukrishnan, and D. Srivastava. Space- and
// time-efficient deterministic algorithms for biased quantiles over data
// streams. In ACM Principles of Database Systems (PODS), 2006.
class SbQuantiles
{
public:
    // Our universe. Inserted items are of this type.
    using U = uint32_t;
    
    // UniverseBits is the number of bits in the universe. Currently it has to
    // be less than 32 (so [0, 2**31-1] is the maximum allowable range of items
    // that can be inserted).
    //
    // epsilon is an accuracy parameter. Smaller epsilon means better accuracy
    // but also a larger data structure. See SbQuantiles::rank for the accuracy
    // guarantee provided.
    //
    // FIXME: SbQuantiles should be a template with UniverseBits as a template
    // parameter.
    SbQuantiles(uint UniverseBits, double epsilon = 0.01);

    // Insert the element 'x'.
    // Precondition: 0 <= x < 2**UniverseBits.
    void insert(U x);
    
    // Compute approximate rank of 'x'.
    // rank(x) is guaranteed to satisfy
    //     abs(rank(x) - true_rank(x)) < epsilon*true_rank(x)
    // where true_rank(x) is the true rank of x.
    int rank(uint64_t x) const;
    
    // Compute approximate 'phi' quantile.
    // Accuracy:
    //    true_rank(x)*(1 − epsilon) <= phi*get_N() <= true_rank(x + 1)*(1 + epsilon)
    // where x = quantile(phi) and true_rank(x) is the true rank of x.
    int quantile(double phi) const;
    
    // Get the number of inserted elements.
    int get_N() const;

    uint64_t get_u_size() const
    {
        return TreeNode::u_size(UniverseBits);
    }

    // The following methods should only be needed during development of
    // SbQuantiles.
    void compress_from_root();
    size_t size() const;
    void dump() const;
    void print_tree() const;
    bool test_invariants() const;
    bool test_accuracy(const ExactQuantiles& e, int seed,
                       uint items_to_test, uint quantiles_to_test) const;
    void set_auto_compress(bool auto_compress);
    void print_params() const;

private:
    using TN = TreeNode;
    struct Node
    {
        Node() : c(0), L(0), w(0)
        { }

        Node(int c_arg, int L_arg) : c(c_arg), L(L_arg), w(0)
        { }

        int c; // Count of this node.
        int L; // Cached L value.
        int w; // Weight of this node. Only used when compressing the tree.
    };
    
    // Number of bits in our universe.
    unsigned UniverseBits;
    
    // The bq-tree.
    using BQ = std::unordered_map<TN, Node>;
    BQ bq;

    double alpha;
    double epsilon;
    
    // Number of items inserted.
    int N;
    
    // Number of items since last compression.
    int num_inserts;

    // If true automatically compress the tree when needed. Will be true during
    // normal operations.
    bool auto_compress;

    int assert_in_u(int x) const;

    void increase_c(TN x, int diff);
    int get_L(TN x) const;
    int get_c(TN x) const;
    int get_c(BQ::const_iterator x) const;
    bool in_bq(TN x) const;
    bool in_bq(BQ::const_iterator x) const;

    int weight(TN v) const;
    void compute_weights();
    int compute_weights_impl(BQ::iterator v_it);
    int get_weight(TN v) const;

    int compute_sum_above(TN x) const;
    void insert_new_leaf(TN leaf);
    void remove_descendants(TN v);
    void compress_tree(TN v, int debt, int L_v, int sum_above);

    int compute_L(TN v) const;
    int compute_A(int x) const;

    void do_auto_compress();
};

#endif
