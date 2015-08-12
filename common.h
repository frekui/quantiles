// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#ifndef common_h
#define common_h

#include <string>
#include <sstream>

#ifndef SBQ_TEST
#  define SBQ_TEST 0
#endif

// The TEST_ASSERT macro is used to verify that certain conditions hold. The
// conditions may be expensive to compute and should only be enabled during
// testing.
//
// Some conditions that always should hold and are cheap to compute use plain
// assert instead.
#if SBQ_TEST
void test_assert_enable();
bool test_assert_is_enabled();
#  define TEST_ASSERT(x) do { if (test_assert_is_enabled()) assert(x); } while(0)
#  ifdef NDEBUG
#    error "SBQ_TEST is true, but NDEBUG is defined"
#  endif
#else
#  define TEST_ASSERT(x)
#endif

template<class C>
std::string to_string(const C& o)
{
    std::stringstream ss;
    ss << o;
    return ss.str();
}

#endif
