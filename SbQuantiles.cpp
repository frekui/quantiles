// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

#include <assert.h>
#include <math.h>

#include "ExactQuantiles.h"
#include "SbQuantiles.h"

using namespace std;

SbQuantiles::SbQuantiles(uint UniverseBits, double epsilon)
{
    assert(UniverseBits <= 31);
    this->UniverseBits = UniverseBits;
    this->epsilon = epsilon;
    auto_compress = true;
    num_inserts = 0;
    N = 0;
    alpha = epsilon/UniverseBits;

    // Set bucket count to a power of two. libcxx will keep it as a power of
    // two.
    // Theory: It will make hashing faster as we replace a modulo computation
    // with a bit mask.
    // Practice: 15% performance improvement as measured by benchmark.cpp.
    bq.rehash(8);
    bq[TN::root()] = Node();
}

void SbQuantiles::print_params() const
{
    cout << "epsilon: " << epsilon
         << " alpha: " << alpha
         << " 1/alpha: " << 1/alpha
         << " |U|: " << TreeNode::u_size(UniverseBits)
         << " UniverseBits: " << UniverseBits
         << " 4/alpha: " << 4/alpha;
}

void SbQuantiles::set_auto_compress(bool auto_compress)
{
    this->auto_compress = auto_compress;
}

int SbQuantiles::assert_in_u(int x) const
{
    assert(0 <= x && x < TreeNode::u_size(UniverseBits));
    return x;
}

int SbQuantiles::get_L(TN x) const
{
    auto it = bq.find(x);
    if (it == bq.end())
        return 0;
    else
        return it->second.L;
}

int SbQuantiles::get_c(TN x) const
{
    return get_c(bq.find(x));
}

int SbQuantiles::get_c(BQ::const_iterator it) const
{
    if (it == bq.end())
        return 0;
    else
        return it->second.c;
}

bool SbQuantiles::in_bq(TN x) const
{
    return in_bq(bq.find(x));
}

bool SbQuantiles::in_bq(BQ::const_iterator it) const
{
    return it != bq.end();
}

void SbQuantiles::increase_c(TN x, int diff)
{
    auto p = bq.insert(make_pair(x, Node(0, 0)));
    auto it = p.first;
    it->second.c += diff;
}

void SbQuantiles::insert_new_leaf(TN leaf)
{
    TEST_ASSERT(leaf.is_leaf(UniverseBits));
    TEST_ASSERT(!in_bq(leaf));

    bq.insert(make_pair(leaf, Node(1, 0)));

    // To maintain the full ancestor property
    // we need to insert ancestors of leaf in bq until we reach
    // an ancestor that is in bq.
    TN x = leaf;
    while (true) {
        x = x.parent(UniverseBits);

        auto p = bq.insert(make_pair(x, Node()));
        // The root is always in bq so we will never go further than that.
        if (!p.second)
            break;
    }
}

// Compress the tree if necessary.
void SbQuantiles::do_auto_compress()
{
    // Condition from Theorem 2.
    if (N < floor(4/alpha))
        return;

    if (N == floor(4/alpha) || num_inserts > log2(epsilon*N)/alpha) {
        // size_t prev_size = size();
        compress_from_root();
        // cout << "Auto-compressing N: " << N << " prev size: " << prev_size << " size: " << size() << endl;
    }
}

void SbQuantiles::insert(U x)
{
    assert_in_u(x);
    N++;
    num_inserts++;
    if (auto_compress)
        do_auto_compress();

    TN leaf = TN::from_universe(UniverseBits, x);

    auto bq_it = bq.find(leaf);
    // If leaf x \in bq, then c_x += 1
    if (bq_it != bq.end()) {
        bq_it->second.c += 1;
        return;
    }

    // Binary search from root to x to find least w \not \in bq.
    TN low = TN::root();
    TN high = leaf;
    // Loop invariant: low \in bq AND up \not \in bq.
    while (low.depth() + 1 < high.depth()) {
        TEST_ASSERT(in_bq(low));
        TEST_ASSERT(!in_bq(high));
        int delta_depth = (high.depth() - low.depth())/2;
        TEST_ASSERT(delta_depth >= 1);
        TEST_ASSERT(high.depth() > delta_depth);
        TN mid = high.parent_k(UniverseBits, delta_depth);
        TEST_ASSERT(mid != low && mid != high);
        if (in_bq(mid)) {
            low = mid;
        } else {
            high = mid;
        }
    }
    // The loop invariant still holds.
    TEST_ASSERT(in_bq(low));
    TEST_ASSERT(!in_bq(high));
    TEST_ASSERT(low.depth() + 1 == high.depth());
    TN w = high;
    TN w_parent = low;
    
    if (get_c(w_parent) + 1 <= floor(alpha*get_L(w_parent))) {
        TEST_ASSERT(in_bq(w_parent));
        increase_c(w_parent, 1);
    } else if (1 <= floor(alpha*get_L(w_parent))) {
        TEST_ASSERT(!in_bq(w));
        bq[w] = Node(1, get_L(w_parent));
    } else {
        insert_new_leaf(leaf);
    }
}

void SbQuantiles::print_tree() const
{
    int depth = 0;
    uint32_t U_size = TreeNode::u_size(UniverseBits);
    int space = 2*U_size;
    cout << depth << ' ';
    TN n = TN::root();
    for (int i = 0; i < 2*U_size - 1; i++) {
        int new_depth = n.depth();
        if (depth != new_depth) {
            TEST_ASSERT(new_depth == depth+1);
            depth++;
            space /= 2;
            cout << '\n' << depth << ' ';
        }
        if (in_bq(n))
            cout << get_c(n);
        else
            cout << '.';
        cout << string(space-1, ' ');
        n = n.bfs_next(UniverseBits);
    }

    cout << endl;
}

int SbQuantiles::weight(TN v) const
{
    int w = get_c(v);
    if (!v.is_leaf(UniverseBits)) {
        if (in_bq(v.left(UniverseBits)))
            w += weight(v.left(UniverseBits));
        if (in_bq(v.right(UniverseBits)))
            w += weight(v.right(UniverseBits));
    }

    return w;
}

void SbQuantiles::compute_weights()
{
    TEST_ASSERT(in_bq(TN::root()));
    compute_weights_impl(bq.find(TN::root()));
}

// Precondition: v_it is a valid iterator in bq.
int SbQuantiles::compute_weights_impl(BQ::iterator v_it)
{
    TEST_ASSERT(v_it != bq.end());
    TN v = v_it->first;
    int w = get_c(v_it);

    if (!v.is_leaf(UniverseBits)) {
        auto v_left_it = bq.find(v.left(UniverseBits));
        if (in_bq(v_left_it))
            w += compute_weights_impl(v_left_it);
        auto v_right_it = bq.find(v.right(UniverseBits));
        if (in_bq(v_right_it))
            w += compute_weights_impl(v_right_it);
    }
    
    v_it->second.w = w;
    return w;
}

void SbQuantiles::remove_descendants(TN v)
{
    if (bq.erase(v) != 0) {
        if (!v.is_leaf(UniverseBits)) {
            remove_descendants(v.left(UniverseBits));
            remove_descendants(v.right(UniverseBits));
        }
    }
}

int SbQuantiles::get_weight(TN v) const
{
    auto it = bq.find(v);
    int w = 0;
    if (in_bq(it))
        w = it->second.w;
    TEST_ASSERT(w == weight(v));
    return w;
}

// Compute the sum \sum_{y ancestor of x and lf(y) = lf(x)} c_y.
// This sum doesn't exist in the paper.
int SbQuantiles::compute_sum_above(TN x) const
{
    if (x.is_root())
        return 0;

    int s = 0;
    TN y = x.parent(UniverseBits);
    while (!y.is_root()) {
        if (x.leftmost_leaf(UniverseBits) != y.leftmost_leaf(UniverseBits))
            break;
        s += get_c(y);
        y = y.parent(UniverseBits);
    }

    return s;
}

static constexpr bool debug_compress_tree = false;

// Implementation of CompressTree in Figure 3 in the paper.
// v:    Start node
// debt: Count borrowed from v and its descendants. Will be payed back by
//       decreasing c_v and/or c_w for some descendants w of v.
// L_v:  L(v)
// sum_above: \sum_{x ancestor of v and lf(v) = lf(x)} c_x
//
// The fourth argument (sum_above) doesn't exist in CompressTree in Figure 3,
// it has been added to fix a problem with the pseudo code in the paper. See
// test.cpp:bug_small for details.
void SbQuantiles::compress_tree(TN v, int debt, int L_v, int sum_above)
{
    if (debug_compress_tree) {
        cout << "v: " << v << " debt: " << debt << " L_v: " << L_v << endl;
        // print_tree();
    }
    TEST_ASSERT(debt >= 0);
    TEST_ASSERT(L_v == compute_L(v));
    TEST_ASSERT(sum_above == compute_sum_above(v));
    auto bq_it = bq.find(v);
    if (bq_it == bq.end()) {
        TEST_ASSERT(debt == 0);
        return;
    }
    Node& node = bq_it->second;
    node.L = L_v;
    int c_prim = node.c;
    if (v.is_leaf(UniverseBits)) {
        node.c = c_prim - debt;
    } else {
        node.c = min(int(floor(alpha*L_v)), get_weight(v) - debt);
        if (node.c >= floor(alpha*L_v)) {
            int new_debt = debt + node.c - c_prim;
            int wl = get_weight(v.left(UniverseBits));

            if (debug_compress_tree)
                cout << string(2*v.depth(), ' ') << "left  ";
            compress_tree(v.left(UniverseBits),
                          min(new_debt, wl),
                          L_v,
                          sum_above + node.c);

            if (debug_compress_tree)
                cout << string(2*v.depth(), ' ') << "right ";
            compress_tree(v.right(UniverseBits),
                          max(new_debt - wl, 0),
                          L_v + wl - min(new_debt, wl) + node.c + sum_above,
                          0);
            // In the paper Line 10 of CompressTree is the following instead
            // compress_tree(v.right(UniverseBits), max(new_debt - wl, 0), L_v + wl + c_prim, 0);
            // However, using the code from the paper seems to sometimes lead to
            // an erronous L(right(v)) value being computed, see
            // test.cpp:bug_small for details.
        } else {
            remove_descendants(v.left(UniverseBits));
            remove_descendants(v.right(UniverseBits));
        }
        if (node.c == 0 &&
            !v.is_root() &&
            !in_bq(v.left(UniverseBits)) &&
            !in_bq(v.right(UniverseBits)))
            bq.erase(v);
    }
}

void SbQuantiles::compress_from_root()
{
    compute_weights();
    compress_tree(TN::root(), 0, 0, 0);
    num_inserts = 0;
}

void SbQuantiles::dump() const
{
    cout << "bq: ";
    for (auto x : bq)
        cout << x.first << ' ';
    cout << endl;
}

// Test the accuracy guarantess for rank and quantiles.
bool SbQuantiles::test_accuracy(const ExactQuantiles& e, int seed,
                              uint items_to_test, uint quantiles_to_test) const
{
    bool ret = true;
    unsigned N = e.get_N();

    std::mt19937 gen;
    gen.seed(seed);
    uint32_t U_size = TreeNode::u_size(UniverseBits);
    std::uniform_int_distribution<> dis(0, U_size-1);
    auto get_element = [&] (int x) {
        if (U_size <= items_to_test)
            return x;
        else
            return dis(gen);
    };
    int test_size = min(U_size, items_to_test);

    if (N != get_N()) {
        cerr << "N = " << N << " != "
             << get_N() << " = get_N()"
             << '\n';
        ret = false;
    }

    // Condition (1), Section 3.1
    // \forall x \in [U]: L(x) - A(x) <= rank(x) <= L(x)
    for (int y = 0; y < test_size; y++) {
        int x = get_element(y);
        int L_x = compute_L(TN::from_universe(UniverseBits, x));
        int A_x = compute_A(x);

        if (L_x - A_x > e.rank(x)) {
            cerr << "(1) failed. x = " << x
                 << "\tL_x - A_x = " << L_x << " - " << A_x << " = " << L_x - A_x
                 << " > " << e.rank(x) << " = rank(x)\n";
            ret = false;
        }
    }

    // Lemma 1, Section 3.1. (Accuracy guarantee)
    for (int y = 0; y < test_size; y++) {
        int x = get_element(y);
        int r_hat = rank(x);
        int r = e.rank(x);

        if (abs(r_hat - r) > epsilon*r) {
            cerr << "Lemma 1 failed. x = " << x
                 << "|r_hat - r| = |" << r_hat << " - " << r
                 << "| = " << abs(r_hat - r)
                 << " > " << epsilon << " = epsilon"
                 << endl;
            ret = false;
        }
    }

    for (double i = 0; i <= quantiles_to_test; i++) {
        double phi = i/double(quantiles_to_test);
        int x_hat = quantile(phi);
        int rank = e.rank(x_hat);
        // Accuracy guarantee for quantiles.
        if (!(rank*(1 - epsilon) <= phi*N)) {
            cerr << "Quantile guarantee failed. "
                 << rank*(1 - epsilon) << " = rank*(1-epsilon) > "
                 << phi*N << " = phi*N\n";
            ret =  false;
        }

        int rank_1 = e.rank(x_hat+1);
        if (!(phi*N <= rank_1*(1 + epsilon))) {
            cerr << "Quantile guarantee failed. "
                 << "phi*N = " << phi*N << " > "
                 << rank_1*(1+epsilon) << " = rank_1*(1+epsilon)\n";
            ret =  false;

        }
    }

    return ret;
}

// Test a number of invariants that the data structure should satisfy. Proofs
// that they hold can be found in the paper.
bool SbQuantiles::test_invariants() const
{
    bool ret = true;

    // In this implementation the root node is always in bq.
    if (!in_bq(TN::root())) {
        cerr << "Root not in bq.";
        ret = false;
    }

    // Condition (2), Section 3.1
    // \forall v \in bq: v \neq lf(v) => c_v \leq floor(\alpha L(v))
    for (auto p : bq) {
        TN v = p.first;
        if (v != v.leftmost_leaf(UniverseBits)) {
            if (get_c(v) > floor(alpha * compute_L(v))) {
                cerr << "(2) failed. v = " << v << " "
                     << get_c(v) << " > floor(" << alpha << "*" << compute_L(v) << ") = "
                     << floor(alpha*compute_L(v)) << endl;
                ret = false;
            }
        }
    }

    // Condition (6), Section 3.2
    // \sum_{v \in bq} c_v = N
    int sum = 0;
    for (auto p : bq) {
        const Node& node = p.second;
        sum += node.c;
    }
    if (sum != N) {
        cerr << "(6) failed. sum = " << sum << " != " << N << " = N\n";
        ret = false;
    }

    // Condition (7), Section 4.1
    // \forall v \in bq : par(v) \in bq and c_par(v) >= floor(alpha L(par(v)))
    // The second part of the condition should hold if num_inserts == 0.
    for (auto p : bq) {
        TN v = p.first;
        if (v.is_root())
            continue;
        TN par = v.parent(UniverseBits);
        if (!in_bq(par)) {
            cerr << "(7) failed. v = " << v
                 << " in bq, but parent(v) = " << v.parent(UniverseBits)
                 << " is not in bq."
                 << endl;
            ret = false;
        }

        if (num_inserts == 0) {
            if (get_c(par) > 0 && get_c(par) < int(floor(alpha*compute_L(par)))) {
                cerr << "(7) failed. par(v) = " << par
                     << " c_par(v) = " << get_c(par)
                     << " < " << floor(alpha*compute_L(par))
                     << endl;
                ret = false;
            }
        }
    }

    // FIXME: Add check for Theorem 1, Section 3.2. (space bound).
    // The theorem needs to be adapted to SBQ, see Theorem 3.

    return ret;
}

// Compute left count, Def. 5
// Tree node -> int.
int SbQuantiles::compute_L(TN v) const
{
    int L_v = 0;
    for (auto p : bq) {
        TN w = p.first;
        if (w.leftmost_leaf(UniverseBits) < v.leftmost_leaf(UniverseBits))
            L_v += get_c(w);
    }
    
    return L_v;
}

// Compute ancestor count, Def. 5
// U -> int
int SbQuantiles::compute_A(int x) const
{
    assert_in_u(x);
    int A_x = 0;
    TN w = TN::from_universe(UniverseBits, x);
    while (!w.is_root()) {
        w = w.parent(UniverseBits);
        A_x += get_c(w);
    }
    A_x += get_c(TN::root());
    return A_x;
}

// We allow x >= u_size as it makes sense and simplifies SbQuantiles::quantile.
int SbQuantiles::rank(uint64_t x) const
{
    if (x >= TreeNode::u_size(UniverseBits))
        return get_N();
    else
        return compute_L(TN::from_universe(UniverseBits, x)) - 0.5*compute_A(x);
}

int SbQuantiles::quantile(double phi) const
{
    assert(0 <= phi && phi <= 1);
    int N = get_N();
    
    // We want to find x \in U so that
    // true_rank(x)*(1 − epsilon) <= phi*N <= true_rank(x + 1)*(1 + epsilon)
    //
    // As true_rank(x)*(1 - epsilon) <= rank(x) <= true_rank(x)*(1 + epsilon)
    // it is sufficient to find x such that rank(x) <= phi*N <= rank(x + 1).

    uint64_t low = 0;
    uint64_t high = TreeNode::u_size(UniverseBits);
    uint64_t x;

    while (true) {
        x = (low + high)/2;
        if (rank(x) > phi*N)
            high = x;
        else if (phi*N > rank(x + 1))
            low = x;
        else
            break;
    }

    return x;
}

size_t SbQuantiles::size() const
{
    return bq.size();
}

int SbQuantiles::get_N() const
{
    return N;
}

