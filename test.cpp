#include <gtest/gtest.h>
#include <fstream>
#include "Functions.h"
#include "IntervalClassesTests.h"
#include "FunctionsTests.h"
#include "ValidationTests.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}