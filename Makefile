CXX = g++
CXXFLAGS = -Wall -g -I/opt/homebrew/Cellar/mpfi/1.5.4/include -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -I/opt/homebrew/Cellar/mpfr/4.2.0-p12/include -I/opt/homebrew/Cellar/googletest/1.14.0/include -std=c++17
LDFLAGS = -L/opt/homebrew/Cellar/mpfi/1.5.4/lib -lmpfi -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -lgmp -L/opt/homebrew/Cellar/mpfr/4.2.0-p12/lib -lmpfr -L/opt/homebrew/Cellar/googletest/1.14.0/lib -lgtest -lgtest_main

all: validation Tests

validation: validation.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

Tests: test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

validation.o: validation.cpp Functions.h IntervalClasses.h
	$(CXX) $(CXXFLAGS) -Wno-c++17-extensions -c $< -o $@

test.o: test.cpp Functions.h IntervalClasses.h IntervalClassesTests.h FunctionsTests.h ValidationTests.h
	$(CXX) $(CXXFLAGS) -Wno-c++17-extensions -c $< -o $@

clean:
	rm -f validation.o validation test.o Tests
