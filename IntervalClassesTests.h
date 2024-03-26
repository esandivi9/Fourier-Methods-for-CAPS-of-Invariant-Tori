#include <gtest/gtest.h>

#include <fstream>

// IntervalTests
// Test if the setInterval function is creating proper intervals around a given value, entered as a string.
TEST(IntervalTest, setInterval) {
    string num = "1.2";
    mpfi_t interval;
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfr_t left, right;
    mpfr_init2(left, prec);
    mpfr_init2(right, prec);
    mpfi_get_left(left, interval);
    mpfi_get_right(right, interval);

    mpfr_mul_si(left, left, 10, MPFR_RNDN);
    mpfr_mul_si(right, right, 10, MPFR_RNDN);

    auto leftSideSmaller = mpfr_cmp_si(left, 12);
    auto rightSideBigger = mpfr_cmp_si(right, 12);

    mpfr_clear(left);
    mpfr_clear(right);
    mpfi_clear(interval);

    EXPECT_TRUE(leftSideSmaller <= 0 && rightSideBigger >= 0);

    //=============================================================

    num = "-2.73";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfr_init2(left, prec);
    mpfr_init2(right, prec);
    mpfi_get_left(left, interval);
    mpfi_get_right(right, interval);

    mpfr_mul_si(left, left, 100, MPFR_RNDN);
    mpfr_mul_si(right, right, 100, MPFR_RNDN);

    leftSideSmaller = mpfr_cmp_si(left, -273);
    rightSideBigger = mpfr_cmp_si(right, -273);

    mpfr_clear(left);
    mpfr_clear(right);
    mpfi_clear(interval);

    EXPECT_TRUE(leftSideSmaller <= 0 && rightSideBigger >= 0);

    //=============================================================

    num = "0.647946395";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfr_init2(left, prec);
    mpfr_init2(right, prec);
    mpfi_get_left(left, interval);
    mpfi_get_right(right, interval);

    mpfr_mul_si(left, left, 1000000000, MPFR_RNDN);
    mpfr_mul_si(right, right, 1000000000, MPFR_RNDN);

    leftSideSmaller = mpfr_cmp_si(left, 647946395);
    rightSideBigger = mpfr_cmp_si(right, 647946395);

    mpfr_clear(left);
    mpfr_clear(right);
    mpfi_clear(interval);

    EXPECT_TRUE(leftSideSmaller <= 0 && rightSideBigger >= 0);

    //=============================================================

    num = "5";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfr_init2(left, prec);
    mpfr_init2(right, prec);
    mpfi_get_left(left, interval);
    mpfi_get_right(right, interval);

    leftSideSmaller = mpfr_cmp_si(left, 5);
    rightSideBigger = mpfr_cmp_si(right, 5);

    mpfr_clear(left);
    mpfr_clear(right);

    EXPECT_TRUE(leftSideSmaller <= 0 && rightSideBigger >= 0);

    //=============================================================

    num = "2.54e2";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfr_init2(left, prec);
    mpfr_init2(right, prec);
    mpfi_get_left(left, interval);
    mpfi_get_right(right, interval);

    leftSideSmaller = mpfr_cmp_si(left, 254);
    rightSideBigger = mpfr_cmp_si(right, 254);

    mpfr_clear(left);
    mpfr_clear(right);

    EXPECT_TRUE(leftSideSmaller <= 0 && rightSideBigger >= 0);

    // Multiples of inverses of powers of 2 should be perfectly represented
    num = "0.125";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfi_t test;
    mpfi_init2(test, prec);
    mpfi_set_d(test, 1. / 8);

    EXPECT_TRUE(mpfi_cmp(test, interval) == 0);

    num = "0.875";
    mpfi_init2(interval, prec);
    setInteval(interval, num);

    mpfi_init2(test, prec);
    mpfi_set_d(test, 7 * 1. / 8);

    EXPECT_TRUE(mpfi_cmp(test, interval) == 0);
}

//==============================================================
// ComplexIntervalTests

TEST(ComplexIntervalTest, AssignationAndConstruction) {
    ComplexInterval a = 1;
    ComplexInterval b(1);
    EXPECT_TRUE(mpfi_cmp(a.real, b.real) == 0);
    EXPECT_TRUE(mpfi_cmp(a.imag, b.imag) == 0);

    a = 1.73;
    ComplexInterval c(1.73);
    EXPECT_TRUE(mpfi_cmp(a.real, c.real) == 0);
    EXPECT_TRUE(mpfi_cmp(a.imag, c.imag) == 0);

    a = "0.56";
    ComplexInterval d("0.56");
    EXPECT_TRUE(mpfi_cmp(a.real, d.real) == 0);
    EXPECT_TRUE(mpfi_cmp(a.imag, d.imag) == 0);
}

TEST(ComplexIntervalTest, ContainedFunction) {
    // Create equal ComplexInterval objects to compare
    ComplexInterval a = "1.0";
    ComplexInterval b = 1.0;

    // Check if one is contained in the other and viceversa
    EXPECT_TRUE(a.isContained(b));
    EXPECT_TRUE(b.isContained(a));

    a = 2;
    b = 3;

    // Check the negation
    EXPECT_TRUE(!a.isContained(b));

    // Create wider intervals (they are big enough to accommodate for floating point error)
    ComplexInterval c;
    ComplexInterval d;

    mpfi_interv_d(c.real, 1.2, 2.4);
    mpfi_interv_d(c.imag, -0.3, 1.7);

    mpfi_interv_d(d.real, 1.3, 2.3);
    mpfi_interv_d(d.imag, -0.2, 1.6);

    EXPECT_TRUE(d.isContained(c));
    EXPECT_FALSE(c.isContained(d));
}

TEST(ComplexIntervalTest, setBounds) {
    ComplexInterval a;
    double realLeft = -1, realRight = -0.125, imagLeft = 3, imagRight = 6;
    setBounds(a, realLeft, realRight, imagLeft, imagRight);
    ComplexInterval expected("-0.5", "5");
    EXPECT_TRUE(expected.isContained(a));
}

TEST(ComplexIntervalTest, AdditionOperator) {
    // Create ComplexInterval objects to add
    ComplexInterval a(1, 1);
    ComplexInterval b(2, 2);

    // Perform the addition
    ComplexInterval result = a + b;

    // Define expected values
    ComplexInterval expected(3, 3);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    ComplexInterval c(1, 2);
    ComplexInterval d(2, 3);

    result = c + d;

    ComplexInterval sum(3, 5);

    EXPECT_TRUE(sum.isContained(result));

    // Check negative decimal numbers
    ComplexInterval e("-1.5", "6.1");
    ComplexInterval f("2.1", "-3.6");

    result = e + f;

    ComplexInterval expSum("0.6", "2.5");

    EXPECT_TRUE(expSum.isContained(result));
}

TEST(ComplexIntervalTest, SubtractionOperator) {
    // Create ComplexInterval objects to add
    ComplexInterval a(1, 0);
    ComplexInterval b(2, 0);

    // Perform the addition
    ComplexInterval result = a - b;

    // Define expected values
    ComplexInterval expected(-1, 0);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    ComplexInterval c(1, 2);
    ComplexInterval d(2, 3);

    result = c - d;

    ComplexInterval sub(-1, -1);

    EXPECT_TRUE(sub.isContained(result));

    // Check negative decimal numbers
    ComplexInterval e("-1.5", "6.1");
    ComplexInterval f("2.1", "-3.6");

    result = e - f;

    ComplexInterval expSub("-3.6", "9.7");

    EXPECT_TRUE(expSub.isContained(result));
}

TEST(ComplexIntervalTest, ProductOperator) {
    // Create ComplexInterval objects to add
    ComplexInterval a(1, 1);
    ComplexInterval b(2, 2);

    // Perform the addition
    ComplexInterval result = a * b;

    // Define expected values
    ComplexInterval expected(0, 4);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    ComplexInterval c(1, 2);
    ComplexInterval d(2, 3);

    result = c * d;

    ComplexInterval prod(-4, 7);

    EXPECT_TRUE(prod.isContained(result));

    // Check negative decimal numbers
    ComplexInterval e("-1.5", "6.1");
    ComplexInterval f("2.1", "-3.6");

    result = e * f;

    ComplexInterval expProd("18.81", "18.21");

    EXPECT_TRUE(expProd.isContained(result));

    //==================================================
    // Multiplication by integer scalar (since double arithmetics are imprecise to test,
    // we will only test integer scalars)
    ComplexInterval interval("2.5", "-3.1");
    double scalar = 2;

    result = interval * scalar;

    ComplexInterval scalarProd("5", "-6.2");

    EXPECT_TRUE(scalarProd.isContained(result));

    // Also multiplication from the left
    result = scalar * interval;
    EXPECT_TRUE(scalarProd.isContained(result));
}

TEST(ComplexIntervalTest, DivisionOperator) {
    // Create ComplexInterval objects to add
    ComplexInterval a(1, 0);
    ComplexInterval b(2, 0);

    // Perform the addition
    ComplexInterval result = a / b;

    // Define expected values
    ComplexInterval expected("0.5", "0");

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    ComplexInterval c(1, -2);
    ComplexInterval d(1, 3);

    result = c / d;

    ComplexInterval div("-0.5", "-0.5");

    EXPECT_TRUE(div.isContained(result));

    // Check other numbers
    ComplexInterval e(1, 2);
    ComplexInterval f(3, 4);

    result = e / f;

    ComplexInterval expDiv("0.44", "0.08");

    EXPECT_TRUE(expDiv.isContained(result));

    // Check decimal numbers
    ComplexInterval g("1.3", "-2.7");
    ComplexInterval h("1.6", "1.6");

    result = g / h;

    ComplexInterval expDivision("-0.4375", "-1.25");

    EXPECT_TRUE(expDivision.isContained(result));

    // Check division by 0
    ComplexInterval zero(0);
    result = a / zero;
    EXPECT_TRUE(mpfi_nan_p(result.real) != 0 || mpfi_nan_p(result.imag) != 0);
}

TEST(ComplexIntervalTest, PowerOperator) {
    // Create ComplexInterval base for exponentiation
    ComplexInterval a(5, 0);

    // Set exponent
    int exp = 1;

    // Perform the addition
    ComplexInterval result = a ^ exp;

    // Define expected values
    ComplexInterval expected(5, 0);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    // Test exponent 0
    exp = 0;
    ComplexInterval expResult(1, 0);
    result = a ^ exp;

    EXPECT_TRUE(expResult.isContained(result));

    // Check other integer values
    ComplexInterval b(-2, 3);
    exp = 3;
    result = b ^ exp;
    ComplexInterval expPower(46, 9);

    EXPECT_TRUE(expPower.isContained(result));

    // Check non-integer values
    ComplexInterval c("1.3", "-3.02");
    exp = 4;
    result = c ^ exp;
    ComplexInterval expExponent("-6.44305984", "116.6870016");

    EXPECT_TRUE(expExponent.isContained(result));
}

TEST(ComplexIntervalTest, OperationAssignationOperator) {
    // Create ComplexInterval base for exponentiation
    ComplexInterval a("1.3", "-2.7");
    ComplexInterval b("1.6", "1.6");

    ComplexInterval copy = a;

    // Addition
    copy += b;
    ComplexInterval expected("2.9", "-1.1");
    EXPECT_TRUE(expected.isContained(copy));
    copy = a;

    // Subtraction
    copy -= b;
    ComplexInterval expResult("-0.3", "-4.3");
    EXPECT_TRUE(expResult.isContained(copy));
    copy = a;

    // Product
    copy *= b;
    ComplexInterval expProd("6.4", "-2.24");
    EXPECT_TRUE(expProd.isContained(copy));
    copy = a;

    // Division
    copy /= b;
    ComplexInterval expDiv("-0.4375", "-1.25");

    EXPECT_TRUE(expDiv.isContained(copy));
}

TEST(ComplexIntervalTest, ComparisonOperator) {
    // Create ComplexInterval base for comparison
    ComplexInterval a("1.3", "-0.7");
    ComplexInterval b("1.3", "-0.7");

    // Equality
    EXPECT_TRUE(a == b);

    // Inequality
    b = 1;
    EXPECT_TRUE(a != b);
}

TEST(ComplexIntervalTest, Sin) {
    // Try the sine of pi/2
    ComplexInterval a = 0;
    mpfi_const_pi(a.real);
    mpfi_div_si(a.real, a.real, 2);

    auto result = sin(a);
    ComplexInterval expResult(1, 0);
    EXPECT_TRUE(expResult.isContained(result));

    ComplexInterval b(1, -1);
    result = sin(b);
    ComplexInterval expSine;
    mpfi_interv_d(expSine.real, 1.29, 1.3);
    mpfi_interv_d(expSine.imag, -0.64, -0.63);

    // Since we made expSine big enough to contain the result, this time the reuslt must be contained in the expected interval
    EXPECT_TRUE(result.isContained(expSine));
}

TEST(ComplexIntervalTest, Cos) {
    // Try the cosine of pi
    ComplexInterval a = 0;
    mpfi_const_pi(a.real);

    auto result = cos(a);
    ComplexInterval expResult(-1, 0);
    EXPECT_TRUE(expResult.isContained(result));

    ComplexInterval b(1, -1);
    result = cos(b);
    ComplexInterval expCosine;
    mpfi_interv_d(expCosine.real, 0.83, 0.84);
    mpfi_interv_d(expCosine.imag, 0.98, 0.99);

    // Since we made expSine big enough to contain the result, this time the reuslt must be contained in the expected interval
    EXPECT_TRUE(result.isContained(expCosine));
}

//==================================================================
// RealIntervalTests

TEST(RealIntervalTest, AssignationAndConstruction) {
    RealInterval a = 1;
    RealInterval b(1);
    EXPECT_TRUE(mpfi_cmp(a.interval, b.interval) == 0);

    a = 1.73;
    RealInterval c(1.73);
    EXPECT_TRUE(mpfi_cmp(a.interval, c.interval) == 0);

    a = "0.56";
    RealInterval d("0.56");
    EXPECT_TRUE(mpfi_cmp(a.interval, d.interval) == 0);
}

TEST(RealIntervalTest, ContainedFunction) {
    // Create equal RealInterval objects to compare
    RealInterval a = "1.0";
    RealInterval b = 1.0;

    // Check if one is contained in the other and viceversa
    EXPECT_TRUE(a.isContained(b));
    EXPECT_TRUE(b.isContained(a));

    a = 2;
    b = 3;

    // Check the negation
    EXPECT_TRUE(!a.isContained(b));

    // Create wider intervals (they are big enough to accommodate for floating point error)
    RealInterval c;
    RealInterval d;

    mpfi_interv_d(c.interval, 1.2, 2.4);
    mpfi_interv_d(d.interval, 1.3, 2.3);

    EXPECT_TRUE(d.isContained(c));
    EXPECT_FALSE(c.isContained(d));
}

TEST(RealIntervalTest, Pi) {
    RealInterval leftB = "3.1415";
    RealInterval rightB = "3.1416";
    EXPECT_TRUE(leftB < PI && PI < rightB);

    leftB = "6.28";
    rightB = "6.29";
    EXPECT_TRUE(leftB < DPI && DPI < rightB);
}

TEST(RealIntervalTest, AdditionOperator) {
    // Create RealInterval objects to add
    RealInterval a(1);
    RealInterval b(2);

    // Perform the addition
    RealInterval result = a + b;

    // Define expected values
    RealInterval expected(3);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    // Check negative decimal numbers
    RealInterval e("-1.5");
    RealInterval f("2.1");

    result = e + f;

    RealInterval expSum("0.6");

    EXPECT_TRUE(expSum.isContained(result));
}

TEST(RealIntervalTest, SubtractionOperator) {
    // Create RealInterval objects to subtract
    RealInterval a(1);
    RealInterval b(2);

    // Perform the addition
    RealInterval result = a - b;

    // Define expected values
    RealInterval expected(-1);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    // Check negative decimal numbers
    RealInterval e("-1.5");
    RealInterval f("2.1");

    result = e - f;

    RealInterval expSub("-3.6");

    EXPECT_TRUE(expSub.isContained(result));
}

TEST(RealIntervalTest, ProductOperator) {
    // Create RealInterval objects to multiply
    RealInterval a(1);
    RealInterval b(2);

    // Perform the addition
    RealInterval result = a * b;

    // Define expected values
    RealInterval expected(2);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    // Check negative decimal numbers
    RealInterval e("-1.5");
    RealInterval f("2.1");

    result = e * f;

    RealInterval expProd("-3.15");

    EXPECT_TRUE(expProd.isContained(result));

    //==================================================
    // Multiplication by integer scalar (since double arithmetics are imprecise to test,
    // we will only test integer scalars)
    RealInterval interval("2.5");
    double scalar = 2;

    result = interval * scalar;

    RealInterval scalarProd("5");

    EXPECT_TRUE(scalarProd.isContained(result));

    // Also multiplication from the left
    result = scalar * interval;
    EXPECT_TRUE(scalarProd.isContained(result));
}

TEST(RealIntervalTest, DivisionOperator) {
    // Create RealInterval objects to divide
    RealInterval a(1);
    RealInterval b(2);

    // Perform the addition
    RealInterval result = a / b;

    // Define expected values
    RealInterval expected("0.5");

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    RealInterval c(1);
    RealInterval d(1);

    result = c / d;

    RealInterval div(1);

    EXPECT_TRUE(div.isContained(result));

    // Check other numbers
    RealInterval e(-7);
    RealInterval f(4);

    result = e / f;

    RealInterval expDiv("-1.75");

    EXPECT_TRUE(expDiv.isContained(result));

    // Check decimal numbers
    RealInterval g("7.3");
    RealInterval h("1.6");

    result = g / h;

    RealInterval expDivision("4.5625");

    EXPECT_TRUE(expDivision.isContained(result));

    // Check division by 0
    RealInterval zero(0);
    result = a / zero;
    EXPECT_TRUE(mpfi_inf_p(result.interval) != 0);
}

TEST(RealIntervalTest, PowerOperator) {
    // Create RealInterval base for exponentiation
    RealInterval a(5);

    // Set exponent
    RealInterval exp = 1;

    // Perform the addition
    RealInterval result = a ^ exp;

    // Define expected values
    RealInterval expected(5);

    // Check if the addition result matches the expected values
    EXPECT_TRUE(expected.isContained(result));

    // Test exponent 0
    exp = 0;
    result = a ^ exp;

    RealInterval expResult(1);

    EXPECT_TRUE(expResult.isContained(result));

    // Check other integer values
    RealInterval b(2);
    exp = 3;
    result = b ^ exp;
    RealInterval expPower(8);

    EXPECT_TRUE(expPower.isContained(result));

    // Check non-integer values
    RealInterval c("1.25");
    exp = -2;
    result = c ^ exp;
    RealInterval expExponent("0.64");

    EXPECT_TRUE(expExponent.isContained(result));
}

TEST(RealIntervalTest, OperationAssignationOperator) {
    RealInterval a("1.3");
    RealInterval b("1.6");

    // Addition
    a += b;
    RealInterval expected("2.9");
    EXPECT_TRUE(expected.isContained(a));
}

TEST(RealIntervalTest, ComparisonOperators) {
    // Create RealInterval base for comparison
    RealInterval a("1.3");
    RealInterval b("1.3");

    // Equality
    EXPECT_TRUE(a == b);
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a >= b);

    // Inequality
    b = 1;
    EXPECT_TRUE(a != b);

    // Order
    EXPECT_TRUE(a > b);
    EXPECT_TRUE(a >= b);
    EXPECT_TRUE(!(b > a));
    EXPECT_TRUE(!(b >= a));
    EXPECT_TRUE(b < a);
    EXPECT_TRUE(b <= a);
    EXPECT_TRUE(!(a < b));
    EXPECT_TRUE(!(a <= b));

    // Overlapping Intervals
    RealInterval c;
    mpfi_interv_d(c.interval, 1.2, 2.4);

    RealInterval d;
    mpfi_interv_d(d.interval, 1, 2);

    EXPECT_TRUE(c > d);
    EXPECT_TRUE(d < c);
    EXPECT_TRUE(!(c < d));
    EXPECT_TRUE(!(d > c));
    EXPECT_TRUE(!(d >= c));
    EXPECT_TRUE(!(c <= d));
}

TEST(RealIntervalTest, SquareRoot) {
    // Create RealInterval base for square root
    RealInterval a(4);
    RealInterval result = sqrt(a);
    RealInterval expResult(2);

    EXPECT_TRUE(expResult.isContained(result));

    a = "24.01";
    result = sqrt(a);
    expResult = "4.9";
    EXPECT_TRUE(expResult.isContained(result));

    // Try approximate roots
    a = 2;
    result = sqrt(a);
    RealInterval lowerB = "1.4";
    RealInterval upperB = "1.5";
    EXPECT_TRUE((lowerB <= result) && (result <= upperB));
}

TEST(RealIntervalTest, AbsoluteValue) {
    RealInterval a(4);
    RealInterval result = abs(a);
    RealInterval expResult(4);

    EXPECT_TRUE(expResult.isContained(result));

    a = "-24.01";
    result = abs(a);
    expResult = "24.01";
    EXPECT_TRUE(expResult.isContained(result));
}

TEST(RealIntervalTest, Sin) {
    // Try the sine of pi/2
    RealInterval a = PI / 2;

    auto result = sin(a);
    RealInterval expResult = 1;
    EXPECT_TRUE(expResult.isContained(result));

    a = 0;
    result = sin(a);
    expResult = 0;
    EXPECT_TRUE(expResult.isContained(result));

    a = PI / 6;
    result = sin(a);
    expResult = "0.5";

    EXPECT_TRUE(expResult.isContained(result));
}

TEST(RealIntervalTest, Cos) {
    // Try the sine of pi/2
    RealInterval a = PI / 2;

    auto result = cos(a);
    RealInterval expResult = 0;
    EXPECT_TRUE(expResult.isContained(result));

    a = 0;
    result = cos(a);
    expResult = 1;
    EXPECT_TRUE(expResult.isContained(result));

    a = PI / 3;
    result = cos(a);
    expResult = "0.5";

    EXPECT_TRUE(expResult.isContained(result));
}

TEST(RealIntervalTest, RealExponential) {
    RealInterval a = 0;
    auto result = exp(a);
    RealInterval expResult = 1;
    EXPECT_TRUE(expResult.isContained(result));

    // Try approximate exp
    a = 1;
    result = exp(a);
    RealInterval lowerB = "2.71";
    RealInterval upperB = "2.72";
    EXPECT_TRUE((lowerB <= result) && (result <= upperB));

    a = -1;
    result = exp(a);
    lowerB = "0.36";
    upperB = "0.37";
    EXPECT_TRUE((lowerB <= result) && (result <= upperB));
}

TEST(RealIntervalTest, Supremum) {
    RealInterval a = 0, b = 1;
    EXPECT_TRUE(abs(b) == sup(a, b));
    EXPECT_TRUE(abs(a) != sup(a, b));

    a = 1;
    EXPECT_TRUE(abs(b) == sup(a, b));
    EXPECT_TRUE(abs(a) == sup(a, b));

    b = -2;
    EXPECT_TRUE(abs(b) == sup(a, b));
    EXPECT_TRUE(abs(a) != sup(a, b));

    a = "3.23";
    b = "-5.43";
    EXPECT_TRUE(abs(b) == sup(a, b));
    EXPECT_TRUE(abs(a) != sup(a, b));
}

//============================================================
// Functions that involve Complex and Real Intervals

TEST(ComplexAndRealIntervalTest, SubtractionOperator) {
    ComplexInterval a("1.5", "-2.7");
    RealInterval b("0.9");

    auto result = a - b;
    ComplexInterval expResult("0.6", "-2.7");
    EXPECT_TRUE(expResult.isContained(result));
}

TEST(ComplexAndRealIntervalTest, ProductOperator) {
    ComplexInterval a("1.5", "-2.7");
    RealInterval b("0.9");

    // RealInterval by ComplexInterval product is thought as a distributive law, so we put the real one in front
    auto result = b * a;
    ComplexInterval expResult("1.35", "-2.43");
    EXPECT_TRUE(expResult.isContained(result));
}

TEST(ComplexAndRealIntervalTest, ComplexModulus) {
    ComplexInterval a(1, 0);

    auto result = mod(a);
    RealInterval expResult(1);
    EXPECT_TRUE(expResult.isContained(result));

    ComplexInterval b(1, -1);
    result = mod(b);
    RealInterval lowerB = "1.4";
    RealInterval upperB = "1.5";
    EXPECT_TRUE((lowerB <= result) && (result <= upperB));
}

TEST(ComplexAndRealIntervalTest, ComplexExponential) {
    RealInterval a = -1 * PI;
    auto result = cExp(a);
    ComplexInterval expResult(-1, 0);
    EXPECT_TRUE(expResult.isContained(result));

    a = PI / 2;
    result = cExp(a);
    ComplexInterval expExponential(0, 1);
    EXPECT_TRUE(expExponential.isContained(result));

    a = PI / 3;
    result = cExp(a);
    ComplexInterval expRes(0.5);
    EXPECT_TRUE(mpfi_is_inside(expRes.real, result.real));

    a = -1 * PI / 6;
    result = cExp(a);
    ComplexInterval expExp(-0.5);
    EXPECT_TRUE(mpfi_is_inside(expExp.imag, result.imag));

    a = 0;
    result = cExp(a);
    ComplexInterval exp(1, 0);
    EXPECT_TRUE(exp.isContained(result));
}

//=========================================================
// Handling vectors and matrices

TEST(VectorComplexIntervalTest, VectorAddition) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 0), 2, ComplexInterval("-0.6", "1.7")};
    vector<ComplexInterval> vec2 = {-3, ComplexInterval("2.6", "0.7"), ComplexInterval(-2, 1)};

    auto result = vec + vec2;

    vector<ComplexInterval> expResult = {ComplexInterval(-2, -3), ComplexInterval("4.6", "2.7"), ComplexInterval("-2.6", "2.7")};

    EXPECT_TRUE(result.size() == expResult.size());

    for (int i = 0; i < vec.size(); i++)
        EXPECT_TRUE(expResult[i].isContained(result[i]));
}

TEST(VectorComplexIntervalTest, VectorSubtraction) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 0), 2, ComplexInterval("-0.6", "1.7")};
    vector<ComplexInterval> vec2 = {-3, ComplexInterval("2.6", "0.7"), ComplexInterval(-2, 1)};

    auto result = vec - vec2;

    vector<ComplexInterval> expResult = {ComplexInterval(4, 3), ComplexInterval("-0.6", "1.3"), ComplexInterval("1.4", "0.7")};

    EXPECT_TRUE(result.size() == expResult.size());

    for (int i = 0; i < vec.size(); i++)
        EXPECT_TRUE(expResult[i].isContained(result[i]));
}

TEST(VectorComplexIntervalTest, VectorEntrywiseProduct) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 0), 2, ComplexInterval("-0.6", "1.7")};
    vector<ComplexInterval> vec2 = {-3, ComplexInterval("2.6", "0.7"), ComplexInterval(-2, 1)};

    auto result = vec * vec2;

    vector<ComplexInterval> expResult = {-3, ComplexInterval("3.8", "6.6"), ComplexInterval("-0.5", "-4")};

    EXPECT_TRUE(result.size() == expResult.size());

    for (int i = 0; i < vec.size(); i++)
        EXPECT_TRUE(expResult[i].isContained(result[i]));
}

//========================================================
// Allocation Tests

TEST(AllocationTest, ComplexMatrixAllocation) {
    unsigned n = 3;
    unsigned m = 4;
    unsigned N = 5;

    ComplexMatrix matrix;
    alloc(matrix, n, m, N);

    // Check if the matrix dimensions match the expected values
    EXPECT_EQ(matrix.size(), n);
    for (unsigned i = 0; i < n; ++i) {
        EXPECT_EQ(matrix[i].size(), m);
        for (unsigned j = 0; j < m; ++j) {
            EXPECT_EQ(matrix[i][j].size(), N);
        }
    }
}

TEST(AllocationTest, ComplexVectorAllocation) {
    unsigned n = 3;
    unsigned N = 5;

    ComplexVector vector;
    alloc(vector, n, N);

    // Check if the vector dimensions match the expected values
    EXPECT_EQ(vector.size(), n);
    for (unsigned i = 0; i < n; ++i) {
        EXPECT_EQ(vector[i].size(), N);
    }
}

TEST(AllocationTest, RealVectorAllocation) {
    unsigned n = 3;
    unsigned N = 5;

    RealVector vector;
    alloc(vector, n, N);

    // Check if the vector dimensions match the expected values
    EXPECT_EQ(vector.size(), n);
    for (unsigned i = 0; i < n; ++i) {
        EXPECT_EQ(vector[i].size(), N);
    }
}

TEST(AllocationTest, ComplexVectorsAllAtOnce) {
    unsigned n = 3;
    unsigned N = 5;
    ComplexVector vector1, vector2, vector3;
    allocAllVectors(n, N, vector1, vector2, vector3);

    // Check if the vector dimensions match the expected values for all vectors
    EXPECT_EQ(vector1.size(), n);
    EXPECT_EQ(vector2.size(), n);
    EXPECT_EQ(vector3.size(), n);
    for (unsigned i = 0; i < n; ++i) {
        EXPECT_EQ(vector1[i].size(), N);
        EXPECT_EQ(vector2[i].size(), N);
        EXPECT_EQ(vector3[i].size(), N);
    }
}

TEST(AllocationTest, ComplexMatricesAllAtOnce) {
    unsigned n = 3;
    unsigned m = 4;
    unsigned N = 5;
    ComplexMatrix matrix1, matrix2, matrix3;
    allocAllMatrices(n, m, N, matrix1, matrix2, matrix3);

    // Check if the matrix dimensions match the expected values for all matrices
    EXPECT_EQ(matrix1.size(), n);
    EXPECT_EQ(matrix2.size(), n);
    EXPECT_EQ(matrix3.size(), n);
    for (unsigned i = 0; i < n; ++i) {
        EXPECT_EQ(matrix1[i].size(), m);
        EXPECT_EQ(matrix2[i].size(), m);
        EXPECT_EQ(matrix3[i].size(), m);
        for (unsigned j = 0; j < m; ++j) {
            EXPECT_EQ(matrix1[i][j].size(), N);
            EXPECT_EQ(matrix2[i][j].size(), N);
            EXPECT_EQ(matrix3[i][j].size(), N);
        }
    }
}

//==================================================
// ComplexMatrix and ComplexVector operations

TEST(ComplexMatrixTest, VectorSubstractionOperator) {
   // Create two ComplexVectors
    ComplexVector vector1, vector2;
    alloc(vector1, 2, 2);
        alloc(vector2, 2, 2);


    vector1[0][0] = ComplexInterval("1.2", "1");
    vector1[1][0] = ComplexInterval("3.1", "-1");
    vector1[0][1] = ComplexInterval("3", "1.5");
    vector1[1][1] = ComplexInterval(2, 3);

    vector2[0][0] = ComplexInterval("1.2", "1");
    vector2[1][0] = ComplexInterval("2.2", "-2");
    vector2[0][1] = ComplexInterval("-1.7", "1.3");
    vector2[1][1] = ComplexInterval(1, 0);

    // Perform matrix multiplication
    ComplexVector result = vector1 - vector2;

    // Define the expected result (manually calculate the expected result)
    ComplexVector expected_result;
    alloc(expected_result, 2, 2);

    expected_result[0][0] = ComplexInterval(0);
    expected_result[1][0] = ComplexInterval("0.9", "1");
    expected_result[0][1] = ComplexInterval("4.7", "0.2");
    expected_result[1][1] = ComplexInterval(1, 3);

    // Check if each element in the result matrix matches the expected result
    for (size_t i = 0; i < expected_result.size(); i++) {
        for (size_t j = 0; j < expected_result[0].size(); j++) {
            // Compare each element
            EXPECT_TRUE(expected_result[i][j].isContained(result[i][j]));
        }
    }
}

TEST(ComplexMatrixTest, MatrixSubstractionOperator) {
    // Create two ComplexMatrices
    ComplexMatrix matrix1;
    alloc(matrix1, 2, 2, 2);

    matrix1[0][0][0] = ComplexInterval(1, -2);
    matrix1[0][1][0] = ComplexInterval(-1, 3);
    matrix1[1][0][0] = ComplexInterval(3, -1);
    matrix1[1][1][0] = ComplexInterval(4, 1);

    matrix1[0][0][1] = ComplexInterval(-1, -2);
    matrix1[0][1][1] = ComplexInterval(-1, 1);
    matrix1[1][0][1] = ComplexInterval(5, -2);
    matrix1[1][1][1] = ComplexInterval(1, 1);

    ComplexMatrix matrix2;
    alloc(matrix2, 2, 2, 2);

    matrix2[0][0][0] = ComplexInterval("1.2", "1");
    matrix2[0][1][0] = ComplexInterval("3", "1.5");
    matrix2[1][0][0] = ComplexInterval("3.1", "-1");
    matrix2[1][1][0] = ComplexInterval(2, 3);

    matrix2[0][0][1] = ComplexInterval("0.2", "1");
    matrix2[0][1][1] = ComplexInterval("0.3", "-1.7");
    matrix2[1][0][1] = ComplexInterval("3.1", "-2");
    matrix2[1][1][1] = ComplexInterval(1, 3);

    // Perform matrix multiplication
    ComplexMatrix result = matrix1 - matrix2;

    // Define the expected result (manually calculate the expected result)
    ComplexMatrix expected_result;
    alloc(expected_result, 2, 2, 2);

    expected_result[0][0][0] = ComplexInterval("-0.2", "-3");
    expected_result[0][1][0] = ComplexInterval("-4", "1.5");
    expected_result[1][0][0] = ComplexInterval("-0.1", "0");
    expected_result[1][1][0] = ComplexInterval("2", "-2");

    expected_result[0][0][1] = ComplexInterval("-1.2", "-3");
    expected_result[0][1][1] = ComplexInterval("-1.3", "2.7");
    expected_result[1][0][1] = ComplexInterval("1.9", "0");
    expected_result[1][1][1] = ComplexInterval("0", "-2");

    // Check if each element in the result matrix matches the expected result
    for (size_t i = 0; i < expected_result.size(); i++) {
        for (size_t j = 0; j < expected_result[0].size(); j++) {
            for (size_t k = 0; k < expected_result[0][0].size(); k++) {
                // Compare each element
                EXPECT_TRUE(expected_result[i][j][k].isContained(result[i][j][k]));
            }
        }
    }
}

TEST(ComplexMatrixTest, MatrixMultiplicationOperator) {
    // Create two ComplexMatrices
    ComplexMatrix matrix1;
    alloc(matrix1, 2, 2, 2);

    matrix1[0][0][0] = ComplexInterval(1, -2);
    matrix1[0][1][0] = ComplexInterval(-1, 3);
    matrix1[1][0][0] = ComplexInterval(3, -1);
    matrix1[1][1][0] = ComplexInterval(4, 1);

    matrix1[0][0][1] = ComplexInterval(-1, -2);
    matrix1[0][1][1] = ComplexInterval(-1, 1);
    matrix1[1][0][1] = ComplexInterval(5, -2);
    matrix1[1][1][1] = ComplexInterval(1, 1);

    ComplexMatrix matrix2;
    alloc(matrix2, 2, 2, 2);

    matrix2[0][0][0] = ComplexInterval("1.2", "1");
    matrix2[0][1][0] = ComplexInterval("3", "1.5");
    matrix2[1][0][0] = ComplexInterval("3.1", "-1");
    matrix2[1][1][0] = ComplexInterval(2, 3);

    matrix2[0][0][1] = ComplexInterval("0.2", "1");
    matrix2[0][1][1] = ComplexInterval("0.3", "-1.7");
    matrix2[1][0][1] = ComplexInterval("3.1", "-2");
    matrix2[1][1][1] = ComplexInterval(1, 3);

    // Perform matrix multiplication
    ComplexMatrix result = matrix1 * matrix2;

    // Define the expected result (manually calculate the expected result)
    ComplexMatrix expected_result;
    alloc(expected_result, 2, 2, 2);

    expected_result[0][0][0] = ComplexInterval("3.1", "8.9");
    expected_result[0][1][0] = ComplexInterval("-5", "-1.5");
    expected_result[1][0][0] = ComplexInterval("18", "0.9");
    expected_result[1][1][0] = ComplexInterval("15.5", "15.5");

    expected_result[0][0][1] = ComplexInterval("0.7", "3.7");
    expected_result[0][1][1] = ComplexInterval("-7.7", "-0.9");
    expected_result[1][0][1] = ComplexInterval("8.1", "5.7");
    expected_result[1][1][1] = ComplexInterval("-3.9", "-5.1");

    // Check if each element in the result matrix matches the expected result
    for (size_t i = 0; i < expected_result.size(); i++) {
        for (size_t j = 0; j < expected_result[0].size(); j++) {
            for (size_t k = 0; k < expected_result[0][0].size(); k++) {
                // Compare each element
                EXPECT_TRUE(expected_result[i][j][k].isContained(result[i][j][k]));
            }
        }
    }
}

TEST(ComplexMatrixTest, MatrixVectorMultiplicationOperator) {
    // Create two ComplexMatrices
    ComplexMatrix matrix;
    alloc(matrix, 2, 2, 2);

    matrix[0][0][0] = ComplexInterval(1, -2);
    matrix[0][1][0] = ComplexInterval(-1, 3);
    matrix[1][0][0] = ComplexInterval(3, -1);
    matrix[1][1][0] = ComplexInterval(4, 1);

    matrix[0][0][1] = ComplexInterval(-1, -2);
    matrix[0][1][1] = ComplexInterval(-1, 1);
    matrix[1][0][1] = ComplexInterval(5, -2);
    matrix[1][1][1] = ComplexInterval(1, 1);

    ComplexVector vector;
    alloc(vector, 2, 2);

    vector[0][0] = ComplexInterval("1.2", "1");
    vector[1][0] = ComplexInterval("3.1", "-1");
    vector[0][1] = ComplexInterval("3", "1.5");
    vector[1][1] = ComplexInterval(2, 3);

    // Perform matrix multiplication
    ComplexVector result = matrix * vector;

    // Define the expected result (manually calculate the expected result)
    ComplexVector expected_result;
    alloc(expected_result, 2, 2);

    expected_result[0][0] = ComplexInterval("3.1", "8.9");
    expected_result[1][0] = ComplexInterval("18", "0.9");

    expected_result[0][1] = ComplexInterval("-5", "-8.5");
    expected_result[1][1] = ComplexInterval("17", "6.5");

    // Check if each element in the result matrix matches the expected result
    for (size_t i = 0; i < expected_result.size(); i++) {
        for (size_t j = 0; j < expected_result[0].size(); j++) {
            // Compare each element
            EXPECT_TRUE(expected_result[i][j].isContained(result[i][j]));
        }
    }
}

TEST(ComplexMatrixTest, MatrixInversion) {
    // Create two ComplexMatrices
    ComplexMatrix matrix;
    alloc(matrix, 2, 2, 2);

    matrix[0][0][0] = ComplexInterval("1.1", "-2");
    matrix[0][1][0] = ComplexInterval("-1.3", "2.1");
    matrix[1][0][0] = ComplexInterval("-3.2", "-0.7");
    matrix[1][1][0] = ComplexInterval(4, 1);

    matrix[0][0][1] = ComplexInterval("4.6", "-0.01");
    matrix[0][1][1] = ComplexInterval("0.87", "1.1");
    matrix[1][0][1] = ComplexInterval(5, -2);
    matrix[1][1][1] = ComplexInterval("9.1", "-3.2");

    // Calculate the inverse
    ComplexMatrix inverseMatrix = inverse(matrix);

    // Perform inverse to the inverse, and check with the original
    ComplexMatrix result = inverse(inverseMatrix);

    // Check if each element in the result matrix matches the expected result
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[0].size(); j++) {
            for (size_t k = 0; k < matrix[0][0].size(); k++) {
                // Compare each element
                EXPECT_TRUE(matrix[i][j][k].isContained(result[i][j][k]));
            }
        }
    }
}