// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#include <assert.h>
#include <unordered_map>

#include "TreeNode.h"

using namespace std;

// Tests for the TreeNode class.

using T = TreeNode;
class TreeNodeTester {
public:
    static int depth(unsigned UniverseBits, T t)
    {
        if (t.x == 0)
            return 0;
        else
            return 1 + depth(UniverseBits, t.parent(UniverseBits));
    }
    
    static T leftmost_leaf(unsigned UniverseBits, T t)
    {
        while (!t.is_leaf(UniverseBits))
            t = t.left(UniverseBits);
        return t;
    }

    static T parent_k(unsigned UniverseBits, T t, int k)
    {
        for (; k > 0; k--) {
            t = t.parent(UniverseBits);
        }

        return t;
    }

    static void test()
    {
        cout << "Running TreeNode tests... " << flush;
        if (!SBQ_TEST)
            abort();

        unsigned UniverseBits = 10;
        T root(UniverseBits, 0);
        TEST_ASSERT(root.left(UniverseBits).x == 1);
        TEST_ASSERT(root.right(UniverseBits).left(UniverseBits).x == 5);
        TEST_ASSERT(T(UniverseBits, 13).parent(UniverseBits).x == 6);
        
        TEST_ASSERT(T::u_size(UniverseBits) == 1024);
        TEST_ASSERT(!T(UniverseBits, 15).is_leaf(UniverseBits));
        TEST_ASSERT(T(UniverseBits, 1024).is_leaf(UniverseBits));
        for (int i = 0; i < T::u_size(UniverseBits)*2-1; i++) {
            T t(UniverseBits, i);
            TEST_ASSERT(T::in_tree(UniverseBits, i));
            TEST_ASSERT(t.depth() == depth(UniverseBits, t));
            TEST_ASSERT(t.leftmost_leaf(UniverseBits) == leftmost_leaf(UniverseBits, t));
            for (int k = 0; k <= t.depth(); k++)
                TEST_ASSERT(t.parent_k(UniverseBits, k) == parent_k(UniverseBits, t, k));
        }
        
        TEST_ASSERT(!T::in_tree(UniverseBits, T::u_size(UniverseBits)*2-1));
        
        unordered_map<T, int> m;
        m[T(UniverseBits, 0)] = 0;
        m[T(UniverseBits, 15)] = 15;
        m[T(UniverseBits, 16)] = 16;
        TEST_ASSERT(m[T(UniverseBits, 0)] == 0);
        TEST_ASSERT(m[T(UniverseBits, 15)] == 15);
        TEST_ASSERT(m[T(UniverseBits, 16)] == 16);
        TEST_ASSERT(m.size() == 3);

        using T3 = TreeNode;
        UniverseBits = 3;
        T3 root3(UniverseBits, 0);
        TEST_ASSERT(root3.leftmost_leaf(UniverseBits).x == root3.left(UniverseBits).left(UniverseBits).left(UniverseBits).x);
        TEST_ASSERT(T3(UniverseBits, 2).leftmost_leaf(UniverseBits).x == T3(UniverseBits, 11).x);
        cout << "ok." << endl;
    }
};

void testTreeNode() {
    TreeNodeTester::test();
}
