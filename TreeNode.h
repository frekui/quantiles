// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#ifndef TreeNode_h
#define TreeNode_h

#include <iostream>
#include <functional>
#include <cstdint>

#include "common.h"

// FIXME: Use a compiler intrinsic instead.
template <class Unsigned>
inline unsigned find_msb(Unsigned mask)
{
    unsigned index = 0;
    while (mask)
    {
        ++index;
        mask >>= 1;
    }
    return --index;
}

class TreeNode;

// A TreeNode represents a node in a complete binary tree using an unsigned
// integer.
//
// The root is represented by 0, the left child of a node 'x' is 2*x+1 and the
// right child is 2*x+2.
//
// FIXME: This class should be a template with UniverseBits as a template
// parameter. UniverseBits is needed in many methods, but as the class isn't a
// template and we don't want to waste space we pass UniverseBits as a parameter
// instead of storing it as a member.
class TreeNode
{
    // Type of universe. Make a better choice later. U should be an unsigned
    // type with at least UniverseBits bits. Each leaf in the tree correponds
    // to one element in U.
    using U = uint32_t;

    // Type of nodes in the tree. The number of bits in T must be strictly
    // greater than the number of bits in U. With the current choice of types
    // U is limited to 31 bits.
    using T = uint32_t;

    explicit TreeNode(unsigned UniverseBits, T a) : x(a)
    {
        assert(in_tree(UniverseBits, x));
    }

public:
    static uint64_t u_size(unsigned UniverseBits)
    {
        return 1ULL << UniverseBits;
    }

    // Convenience.
    using TN = TreeNode;
    
    // Convert an element in U to a TreeNode (a leaf).
    static TN from_universe(unsigned UniverseBits, U u)
    {
        assert(u < (U(1) << UniverseBits));
        return TN(UniverseBits, u + u_size(UniverseBits) - 1);
    }
    
    static TN root()
    {
        return TN(1, 0);
    }

    TN parent(unsigned UniverseBits) const
    {
        return TN(UniverseBits, (x-1)/2);
    }

    // Equivalent to parent(parent(...)) iterated k times.
    TN parent_k(unsigned UniverseBits, int k) const
    {
        return TN(UniverseBits, (x - (1 << k) + 1) >> k);
    }

    static bool in_tree(unsigned UniverseBits, T x)
    {
        return x <= u_size(UniverseBits)*2 - 2;
    }
    
    bool is_root() const
    {
        return x == 0;
    }
    
    bool is_leaf(unsigned UniverseBits) const
    {
        return u_size(UniverseBits) - 2 < x;
    }

    int depth() const
    {
        if (x == 0)
            return 0;
        
        return find_msb(x + 1);
    }

    TN left(unsigned UniverseBits) const
    {
        return TN(UniverseBits, 2*x + 1);
    }

    TN right(unsigned UniverseBits) const
    {
        return TN(UniverseBits, 2*x + 2);
    }
    
    // Named lf(x) in the paper, section 3.
    TN leftmost_leaf(unsigned UniverseBits) const
    {
        unsigned tree_height = UniverseBits;
        return left_k(UniverseBits, tree_height - depth());
    }

    // Equivalent to left(left(...)) iterated k times.
    TN left_k(unsigned UniverseBits, int k) const
    {
        return TN(UniverseBits, (x << k) + ((1 << k) - 1));
    }

    // Return next node in the tree in a BFS.
    // We don't allow invalid TreeNodes so at the final node we wrap around to
    // the root instead of e.g. returning a one-past-the-end TreeNode.
    TN bfs_next(unsigned UniverseBits) const
    {
        if (in_tree(UniverseBits, x + 1))
            return TN(UniverseBits, x + 1);
        else
            return TN::root();
    }

    bool operator==(const TN &rhs) const
    {
        return x == rhs.x;
    }

    bool operator!=(const TN &rhs) const
    {
        return !(*this == rhs);
    }
    
    bool operator<(const TN &rhs) const
    {
        return x < rhs.x;
    }

    friend class TreeNodeTester;
    friend class std::hash<TN>;
    friend std::ostream& operator<< (std::ostream&, const TN&);

private:
    T x;
};

inline std::ostream& operator<< (std::ostream& o, const TreeNode& tn)
{
    return o << tn.x;
}

namespace std {
    template<> class hash<TreeNode>
    {
        std::hash<typename TreeNode::T> h;

    public:
        size_t operator() (const TreeNode& t) const
        {
            return h(t.x);
        }
    };
}
void testTreeNode();

#endif
