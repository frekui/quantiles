// Copyright (c) 2015 Fredrik Kuivinen, frekui@gmail.com.

#include "common.h"

static bool enabled = false;
void test_assert_enable()
{
    enabled = true;
}

bool test_assert_is_enabled()
{
    return enabled;
}
