all: example benchmark test test-big test-cov test-big-cov

CXXFLAGS=-Wall -Wextra -std=gnu++11 -MMD -MP -Werror $(EXTRA_CXXFLAGS)
OPTFLAGS=-O2 -flto
COVFLAGS=-g -O0 --coverage
LINK.o=$(CXX) $(OPTFLAGS)
COMPILE.cpp=$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c
OBJS=ExactQuantiles.o SbQuantiles.o TreeNode.o common.o
OBJS_TEST=$(patsubst %.o, %-test.o, $(OBJS))
OBJS_COV=$(patsubst %.o, %-cov.o, $(OBJS))

ifeq ($(shell $(CXX) --version | grep -q clang && echo y),y)
    $(info CXX identified as clang)
    EXTRA_CXXFLAGS=-Wno-mismatched-tags -Wno-sign-compare
    ifeq ($(shell echo 'int main(void) { return 0; }' | $(CXX) $(CXXFLAGS) $(OPTFLAGS) -x c++ - -o minimal || echo n),n)
        $(error Failed to compile and link minimal C++ program. Clang with -flto not usable?)
    endif
else ifeq ($(shell $(CXX) --version | grep -q g++ && echo y),y)
    $(info CXX identified as gcc)
    EXTRA_CXXFLAGS=-Wno-sign-compare
else
    $(info CXX not identified, things may or may not work.)
endif

-include $(wildcard *.d)

example: $(OBJS) example.o

benchmark: $(OBJS) benchmark.o

# SBG_DEBUG=1 enables lots of fairly expensive assertions that are
# used for testing.
%-test.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -DSBQ_TEST=1 -c $< -o $@

%-cov.o: %.cpp
	$(CXX) $(CXXFLAGS) $(COVFLAGS) -DSBQ_TEST=1 -c $< -o $@

test: $(OBJS_TEST) test-test.o
	$(LINK.o) $(filter %.o, $^) -o $@

test-cov: $(OBJS_COV) test-cov.o
	$(CXX) $(COVFLAGS) $(filter %.o, $^) -o $@

test-big: $(OBJS_TEST) test-big-test.o
	$(LINK.o) $(OPTFLAGS) $(filter %.o, $^) -o $@

test-big-cov: $(OBJS_COV) test-big-cov.o
	$(LINK.o) $(COVFLAGS) $(filter %.o, $^) -o $@

# Run all tests.
.PHONY: run-tests
run-tests: test test-big
	./test
	./test-big

# Run all tests and generate coverage report.
.PHONY: run-coverage
# FIXME: Use gcovr instead?
run-coverage: test-cov test-big-cov
	./test-cov
	./test-big-cov
	gcov TreeNode-cov.h SbQuantiles-cov.cpp SbQuantiles-cov.h test-cov.cpp test-big-cov.cpp
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory covout
	$(info Coverage report written to covout/index.html)
	-open covout/index.html

.PHONY: clean
clean:
	rm -f *~ *.o *.d *.gcda *.gcno
	rm -f test test-big test-cov test-big-cov benchmark example coverage minimal coverage.info
	rm -rf covout
