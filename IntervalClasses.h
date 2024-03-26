#include <mpfi.h>
#include <mpfi_io.h>

#include <complex>
#include <iostream>
#include <numeric>

using namespace std;
mp_prec_t prec = 100;

/*  Since float point representation is not precise, the intervals created around
    such values may not include the values themselves. So we create the interval
    by receiving the desired value in a string, then moving the decimal point to the right,
    creating the interval with such an integer, and dividing it by the interval of as many powers of 10
    as many decimal points we had.
*/
void setInteval(mpfi_t &interval, string num) {
    size_t exponentPos = num.find("e");
    size_t decimalPos = num.find(".");
    int count = 0, extendedNum, power;

    // Check if a decimal point was found
    if (decimalPos != string::npos) {
        // Check if it's written in scientific notation
        if (exponentPos != string::npos) {
            int decimals = num.substr(decimalPos + 1, exponentPos - (decimalPos + 1)).length();  // Check for values between . and e
            auto exp = num.substr(exponentPos + 1);                                              // Extract the exponent part
            count = atoi(exp.c_str());
            count -= decimals;
        } else {
            // Calculate the number of decimal places
            count = -(num.length() - decimalPos - 1);
        }

        num.erase(remove(num.begin(), num.end(), '.'), num.end());
    } else if (exponentPos != string::npos) {
        // No decimal point but scientific notation found
        auto exp = num.substr(exponentPos + 1);  // Extract the exponent part
        count = atoi(exp.c_str());
    }

    extendedNum = stoi(num);
    power = pow(10, abs(count));

    mpfi_t extendedInterval;
    mpfi_init2(extendedInterval, prec);
    mpfi_set_d(extendedInterval, extendedNum);

    mpfi_t powerInterval;
    mpfi_init2(powerInterval, prec);
    mpfi_set_d(powerInterval, power);

    if (count < 0)
        mpfi_div(interval, extendedInterval, powerInterval);
    else if (count > 0)
        mpfi_mul(interval, extendedInterval, powerInterval);
    else
        mpfi_set(interval, extendedInterval);
}

// Custom class for representing complex intervals
class ComplexInterval {
   public:
    mpfi_t real;  // Real part of the interval
    mpfi_t imag;  // Imaginary part of the interval

    // Constructor
    ComplexInterval() {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
    }

    // Constructor that accepts a string
    ComplexInterval(const string &x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        setInteval(real, x);
        setInteval(imag, x);
    }

    // Constructor that accepts a const char* (useful for assignation operator)
    ComplexInterval(const char *x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        setInteval(real, x);
        setInteval(imag, x);
    }

    ComplexInterval(const string &x, const string &y) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        setInteval(real, x);
        setInteval(imag, y);
    }

    ComplexInterval(const complex<double> &x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        mpfi_set_d(real, std::real(x));
        mpfi_set_d(imag, std::imag(x));
    }

    ComplexInterval(const double &x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        mpfi_set_d(real, x);
        mpfi_set_d(imag, x);
    }

    ComplexInterval(const int &x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        mpfi_set_si(real, x);
        mpfi_set_si(imag, x);
    }

    ComplexInterval(const double &x, const double &y) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        mpfi_set_d(real, x);
        mpfi_set_d(imag, y);
    }

    ComplexInterval(const ComplexInterval &x) {
        mpfi_init2(real, prec);
        mpfi_init2(imag, prec);
        mpfi_set(real, x.real);
        mpfi_set(imag, x.imag);
    }

    // Destructor
    ~ComplexInterval() {
        mpfi_clear(real);
        mpfi_clear(imag);
    }

    ComplexInterval &operator=(const ComplexInterval &x) {
        if (this != &x) {
            // Clear existing data
            mpfi_clear(real);
            mpfi_clear(imag);

            // Perform deep copy
            mpfi_init2(real, prec);
            mpfi_init2(imag, prec);
            mpfi_set(real, x.real);
            mpfi_set(imag, x.imag);
        }
        return *this;
    }

    bool isContained(const ComplexInterval &other) const {
        // Check if the real and imaginary parts of 'this' are contained within 'other'
        return (mpfi_is_inside(real, other.real)) && (mpfi_is_inside(imag, other.imag));
    }

    // Overload the + operator for adding two complex intervals
    ComplexInterval operator+(const ComplexInterval &other) const {
        ComplexInterval result;
        mpfi_add(result.real, real, other.real);
        mpfi_add(result.imag, imag, other.imag);
        return result;
    }

    // Overload the - operator for subtracting two complex intervals
    ComplexInterval operator-(const ComplexInterval &other) const {
        ComplexInterval result;
        mpfi_sub(result.real, real, other.real);
        mpfi_sub(result.imag, imag, other.imag);
        return result;
    }

    // Overload the * operator for multiplying two complex intervals
    ComplexInterval operator*(const ComplexInterval &other) const {
        ComplexInterval result;
        mpfi_t temp1, temp2, temp3;
        mpfi_inits(temp1, temp2, temp3, NULL);

        // Real part calculation
        mpfi_mul(temp1, real, other.real);
        mpfi_mul(temp2, imag, other.imag);
        mpfi_sub(result.real, temp1, temp2);

        // Imaginary part calculation
        mpfi_mul(temp3, real, other.imag);
        mpfi_mul(temp2, imag, other.real);
        mpfi_add(result.imag, temp3, temp2);

        mpfi_clears(temp1, temp2, temp3, NULL);
        return result;
    }

    // Overload the * operator for multiplying a complex interval with a scalar
    ComplexInterval operator*(double a) const {
        ComplexInterval z;
        mpfi_mul_d(z.real, real, a);
        mpfi_mul_d(z.imag, imag, a);
        return z;
    }

    ComplexInterval operator/(const ComplexInterval &y) const {
        ComplexInterval result;
        mpfi_t numeratorReal, numeratorImag, denominator, temp;
        mpfi_inits2(prec, numeratorReal, numeratorImag, denominator, temp, NULL);

        // Numerator (real part)
        mpfi_mul(numeratorReal, real, y.real);
        mpfi_mul(numeratorImag, imag, y.imag);
        mpfi_add(result.real, numeratorReal, numeratorImag);

        // Denominator
        mpfi_mul(temp, y.real, y.real);
        mpfi_mul(denominator, y.imag, y.imag);
        mpfi_add(denominator, temp, denominator);

        // Division of real part
        mpfi_div(result.real, result.real, denominator);

        // Numerator (imaginary part)
        mpfi_mul(numeratorReal, real, y.imag);
        mpfi_mul(numeratorImag, imag, y.real);
        mpfi_sub(result.imag, numeratorImag, numeratorReal);

        // Division of imaginary part
        mpfi_div(result.imag, result.imag, denominator);

        mpfi_clears(numeratorReal, numeratorImag, denominator, temp, NULL);
        return result;
    }

    // Overload the ^ operator for exponentiation of a complex interval to an integer power
    ComplexInterval operator^(int n) const {
        ComplexInterval result;
        mpfi_t temp1, temp2, temp3;
        mpfi_inits2(prec, temp1, temp2, temp3, NULL);

        // Initialize result with identity element
        mpfi_set_d(result.real, 1.0);
        mpfi_set_d(result.imag, 0.0);

        // Exponentiation
        for (int i = 0; i < n; ++i) {
            mpfi_mul(temp1, result.real, real);
            mpfi_mul(temp2, result.imag, imag);
            mpfi_sub(temp3, temp1, temp2);

            mpfi_mul(temp1, result.real, imag);
            mpfi_mul(temp2, result.imag, real);
            mpfi_add(result.imag, temp1, temp2);

            mpfi_set(result.real, temp3);
        }

        mpfi_clears(temp1, temp2, temp3, NULL);
        return result;
    }

    // Overload of the += operator
    ComplexInterval &operator+=(const ComplexInterval &other) {
        *this = *this + other;
        return *this;
    }

    // Overload of the -= operator
    ComplexInterval &operator-=(const ComplexInterval &other) {
        *this = *this - other;
        return *this;
    }

    // Overload of the *= operator
    ComplexInterval &operator*=(const ComplexInterval &other) {
        *this = *this * other;
        return *this;
    }

    // Overload of the /= operator
    ComplexInterval &operator/=(const ComplexInterval &other) {
        *this = *this / other;
        return *this;
    }

    // Overload comparison operator ==
    bool operator==(const ComplexInterval &other) const {
        return (mpfi_cmp(real, other.real) == 0) && (mpfi_cmp(imag, other.imag) == 0);
    }

    // Overload comparison operator !=
    bool operator!=(const ComplexInterval &other) const {
        return (mpfi_cmp(real, other.real) != 0) || (mpfi_cmp(imag, other.imag) != 0);
    }
};

// Set the bounds of the intervals for each real and imaginary part
void setBounds(ComplexInterval &x, const double &realLeft,
               const double &realRight, const double &imagLeft,
               const double &imagRight) {
    mpfi_interv_d(x.real, realLeft, realRight);
    mpfi_interv_d(x.imag, imagLeft, imagRight);
}

// Overload the * operator for multiplying a complex interval with a scalar on the left side
ComplexInterval operator*(double a, const ComplexInterval &z) {
    return z * a;
}

// Overload the + operator for vectors of ComplexInterval
vector<ComplexInterval> operator+(const vector<ComplexInterval> &a, const vector<ComplexInterval> &b) {
    vector<ComplexInterval> result;
    result.reserve(a.size());  // Reserve space for the result

    for (size_t i = 0; i < a.size(); i++) {
        result.push_back(a[i] + b[i]);
    }

    return result;
}

// Overload the - operator for vectors of ComplexInterval
vector<ComplexInterval> operator-(const vector<ComplexInterval> &a, const vector<ComplexInterval> &b) {
    vector<ComplexInterval> result;
    result.reserve(a.size());  // Reserve space for the result

    for (size_t i = 0; i < a.size(); i++) {
        result.push_back(a[i] - b[i]);
    }

    return result;
}

// Entry-wise product of vectors of complex intervals
vector<ComplexInterval> operator*(const vector<ComplexInterval> &a, const vector<ComplexInterval> &b) {
    vector<ComplexInterval> result;
    result.reserve(a.size());  // Reserve space for the result

    for (size_t i = 0; i < a.size(); i++) {
        result.push_back(a[i] * b[i]);
    }

    return result;
}

void print(const ComplexInterval &a) {
    mpfi_out_str(stdout, 10, 0, a.real);
    cout << " + i";
    mpfi_out_str(stdout, 10, 0, a.imag);
    cout << endl;
}

class RealInterval {
   public:
    mpfi_t interval;  // The interval

    // Constructor
    RealInterval() {
        mpfi_init2(interval, prec);
    }

    // Constructor that accepts a string
    RealInterval(const string &x) {
        mpfi_init2(interval, prec);
        setInteval(interval, x);
    }

    // Constructor that accepts a const char* (useful for assignation operator)
    RealInterval(const char *x) {
        mpfi_init2(interval, prec);
        setInteval(interval, x);
    }

    RealInterval(const double &x) {
        mpfi_init2(interval, prec);
        mpfi_set_d(interval, x);
    }

    RealInterval(const int &x) {
        mpfi_init2(interval, prec);
        mpfi_set_si(interval, x);
    }

    RealInterval(const RealInterval &x) {
        mpfi_init2(interval, prec);
        mpfi_set(interval, x.interval);
    }

    // Destructor
    ~RealInterval() {
        mpfi_clear(interval);
    }

    // Overload the = operator
    RealInterval &operator=(const RealInterval &x) {
        if (this != &x) {
            // Clear existing data
            mpfi_clear(interval);

            // Perform deep copy
            mpfi_init2(interval, prec);
            mpfi_set(interval, x.interval);
        }
        return *this;
    }

    bool isContained(const RealInterval &other) const {
        // Check if 'this' is contained within 'other'
        return mpfi_is_inside(interval, other.interval);
    }

    // Overload the + operator for adding two real intervals
    RealInterval operator+(const RealInterval &other) const {
        RealInterval result;
        mpfi_add(result.interval, interval, other.interval);
        return result;
    }

    // Overload the - operator for subtracting two real intervals
    RealInterval operator-(const RealInterval &other) const {
        RealInterval result;
        mpfi_sub(result.interval, interval, other.interval);
        return result;
    }

    // Overload the * operator for multiplying two real intervals
    RealInterval operator*(const RealInterval &other) const {
        RealInterval result;
        mpfi_mul(result.interval, interval, other.interval);
        return result;
    }

    // Overload the * operator for multiplying a real interval by a scalar
    RealInterval operator*(const double &other) const {
        RealInterval result;
        mpfi_mul_d(result.interval, interval, other);
        return result;
    }

    // Overload the / operator for dividing two real intervals
    RealInterval operator/(const RealInterval &other) const {
        RealInterval result;
        mpfi_div(result.interval, interval, other.interval);
        return result;
    }

    // Overload the ^ operator for exponentiating a real interval to a real interval power
    RealInterval operator^(const RealInterval &exponent) const {
        RealInterval result;
        if (mpfi_cmp_si(interval, 0) < 0) {
            cout << "Real powers of real negative numbers not allowed. Closing program." << endl;
            exit(0);
        }

        mpfi_log(result.interval, interval);                            // Compute log(|interval|)
        mpfi_mul(result.interval, result.interval, exponent.interval);  // Multiply by the exponent
        mpfi_exp(result.interval, result.interval);                     // Compute exp(result) to get the final result

        return result;
    }

    // Overload of the += operator
    RealInterval &operator+=(const RealInterval &other) {
        *this = *this + other;
        return *this;
    }

    bool operator==(const RealInterval &other) const {
        return mpfi_cmp(interval, other.interval) == 0;
    }

    bool operator!=(const RealInterval &other) const {
        return mpfi_cmp(interval, other.interval) != 0;
    }

    // Comparison operators
    bool operator>(const RealInterval &other) const {
        mpfr_t rightBoundThis, rightBoundOther;
        mpfr_init2(rightBoundThis, prec);
        mpfr_init2(rightBoundOther, prec);
        mpfi_get_right(rightBoundThis, interval);
        mpfi_get_right(rightBoundOther, other.interval);
        int cmpResult = mpfr_cmp(rightBoundThis, rightBoundOther);
        return cmpResult > 0;
    }

    bool operator<(const RealInterval &other) const {
        return other > *this;
    }

    bool operator>=(const RealInterval &other) const {
        return !(*this < other);
    }

    bool operator<=(const RealInterval &other) const {
        return !(other < *this);
    }
};

// Overload the * operator for multiplying a real interval with a scalar on the left side
RealInterval operator*(double a, const RealInterval &x) {
    return x * a;
}

// Overload the - operator for substracting a real interval from a complex interval
ComplexInterval operator-(const ComplexInterval &x, const RealInterval &y) {
    ComplexInterval result = x;
    mpfi_sub(result.real, x.real, y.interval);
    return result;
}

// Multiply a real interval by a complex interval
ComplexInterval operator*(const RealInterval &a, const ComplexInterval &x) {
    ComplexInterval z;
    mpfi_mul(z.real, a.interval, x.real);
    mpfi_mul(z.imag, a.interval, x.imag);
    return z;
}

void print(const RealInterval &a) {
    mpfi_out_str(stdout, 10, 0, a.interval);
    cout << endl;
}

// Positive square root of a real interval
RealInterval sqrt(const RealInterval &x) {
    RealInterval sqrt;
    mpfi_sqrt(sqrt.interval, x.interval);
    return sqrt;
}

// Compute the absolute value of a real interval
RealInterval abs(const RealInterval &x) {
    RealInterval abs;
    mpfi_abs(abs.interval, x.interval);
    return abs;
}

// Calculate the modulus of a complex interval
RealInterval mod(const ComplexInterval &z) {
    RealInterval mod;
    mpfi_hypot(mod.interval, z.real, z.imag);
    return mod;
}

// Compute complex sine of a complex interval
ComplexInterval sin(const ComplexInterval &x) {
    ComplexInterval z = 0;
    RealInterval a, b;

    mpfi_sin(a.interval, x.real);
    mpfi_cosh(b.interval, x.imag);
    mpfi_mul(z.real, a.interval, b.interval);

    mpfi_cos(a.interval, x.real);
    mpfi_sinh(b.interval, x.imag);
    mpfi_mul(z.imag, a.interval, b.interval);

    return z;
}

// Compute complex cosine of a complex interval
ComplexInterval cos(const ComplexInterval &x) {
    ComplexInterval z = 0;
    RealInterval a, b, mOne = -1.0;

    mpfi_cos(a.interval, x.real);
    mpfi_cosh(b.interval, x.imag);
    mpfi_mul(z.real, a.interval, b.interval);

    mpfi_sin(a.interval, x.real);
    mpfi_sinh(b.interval, x.imag);
    mpfi_mul(b.interval, mOne.interval, b.interval);
    mpfi_mul(z.imag, a.interval, b.interval);

    return z;
}

// Calculate the sine of a real interval
RealInterval sin(const RealInterval &x) {
    RealInterval sine;
    mpfi_sin(sine.interval, x.interval);
    return sine;
}

// Calculate the cosine of a real interval
RealInterval cos(const RealInterval &x) {
    RealInterval cosine;
    mpfi_cos(cosine.interval, x.interval);
    return cosine;
}

// Calculate the complex exponential with real interval exponent, e^ix
ComplexInterval cExp(const RealInterval &x) {
    ComplexInterval exp;
    mpfi_cos(exp.real, x.interval);
    mpfi_sin(exp.imag, x.interval);
    return exp;
}

// Calculate the exponential of a real interval, e^x
RealInterval exp(const RealInterval &x) {
    RealInterval exp;
    mpfi_exp(exp.interval, x.interval);
    return exp;
}

// Find the supremum of two real intervals
RealInterval sup(const RealInterval &x, const RealInterval &y) {
    RealInterval sup1 = abs(x);
    RealInterval sup2 = abs(y);

    return (sup1 <= sup2) ? sup2 : sup1;
}

// Define Pi globally using the MPFI pi builder
RealInterval initializePiInterval() {
    mpfi_t pi_interval;
    mpfi_init2(pi_interval, prec);
    mpfi_const_pi(pi_interval);
    RealInterval pi;
    mpfi_set(pi.interval, pi_interval);
    return pi;
}

RealInterval PI = initializePiInterval();
RealInterval DPI = 2 * PI;
RealInterval one = RealInterval(1);
RealInterval two = RealInterval(2);

/* Although using three nested vectors suggests a 3 dimensional matrix,
we use the last dimension as a coefficient array. That is, a
ComplexMatrix will be a matrix in which each entry we have an array
of intervals, and a ComplexVector will be a vector in which each
entry is an array of intervals. The same goes for RealMatrix and others. */

using ComplexMatrix = vector<vector<vector<ComplexInterval>>>;
using ComplexVector = vector<vector<ComplexInterval>>;
using RealVector = vector<vector<RealInterval>>;

// Allocate a matrix of complex intervals
void alloc(ComplexMatrix &M, unsigned n, unsigned m, unsigned N) {
    M.resize(n, ComplexVector(m, vector<ComplexInterval>(N)));
}

// Allocate a set of ComplexMatrices at once
template <typename... Matrices>
void allocAllMatrices(unsigned n, unsigned m, unsigned N, Matrices &...matrices) {
    (alloc(matrices, n, m, N), ...);
}

// Allocate a vector of complex intervals
void alloc(ComplexVector &v, unsigned n, unsigned N) {
    v.resize(n, vector<ComplexInterval>(N));
}

// Allocate a set of ComplexVectors at once
template <typename... Vectors>
void allocAllVectors(unsigned n, unsigned N, Vectors &...vectors) {
    (alloc(vectors, n, N), ...);
}

// Allocate a vector of real intervals
void alloc(RealVector &v, unsigned n, unsigned N) {
    v.resize(n, vector<RealInterval>(N));
}

// Overload the - operator for ComplexVectors
ComplexVector operator-(const ComplexVector &a, const ComplexVector &b) {
    ComplexVector result;
    int m = a.size(), N = a[0].size();

    alloc(result, m, N);

    for (int i = 0; i < m; i++)
        for (int k = 0; k < N; k++)
            result[i][k] = a[i][k] - b[i][k];

    return result;
}

// Overload the - operator for ComplexMatrices
ComplexMatrix operator-(const ComplexMatrix &a, const ComplexMatrix &b) {
    ComplexMatrix result;
    int m = a.size(), N = a[0][0].size();

    alloc(result, m, m, N);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            for (int k = 0; k < N; k++)
                result[i][j][k] = a[i][j][k] - b[i][j][k];

    return result;
}

// Overload * operator for matrix multiplication
ComplexMatrix operator*(const ComplexMatrix &matrix1, const ComplexMatrix &matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();
    int N = matrix1[0][0].size();

    ComplexMatrix result;
    alloc(result, rows1, cols2, N);

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < rows1; i++) {
            for (int j = 0; j < cols2; j++) {
                result[i][j][k] = 0;
                for (int l = 0; l < cols1; l++) {
                    result[i][j][k] += matrix1[i][l][k] * matrix2[l][j][k];
                }
            }
        }
    }

    return result;
}

// Overload * operator for matrix-vector multiplication
ComplexVector operator*(const ComplexMatrix &matrix, const ComplexVector &vector) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    int N = vector[0].size();

    ComplexVector result;
    alloc(result, rows, N);

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < rows; i++) {
            result[i][k] = 0;
            for (int j = 0; j < cols; j++) {
                result[i][k] += matrix[i][j][k] * vector[j][k];
            }
        }
    }

    return result;
}

// Function to swap two rows in a matrix
void swapRows(ComplexMatrix &matrix, int row1, int row2) {
    ComplexVector temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
}

// Function to find the inverse of a matrix of complex intervals
ComplexMatrix inverse(const ComplexMatrix &matrix) {
    ComplexInterval zero(0);
    // Get the dimensions of the matrix
    int n = matrix.size();
    int N = matrix[0][0].size();

    // Create an augmented matrix [matrix | identity]
    ComplexMatrix augmentedMatrix;
    alloc(augmentedMatrix, n, 2 * n, N);
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 2 * n; ++j) {
                augmentedMatrix[i][j][k] = 0;
                if (j < n)
                    augmentedMatrix[i][j][k] = matrix[i][j][k];
            }
            augmentedMatrix[i][i + n][k] = ComplexInterval(1, 0);
        }
    }

    for (int k = 0; k < N; k++) {
        // Apply Gauss-Jordan elimination
        for (int i = 0; i < n; ++i) {
            // Swap rows if the current row has a zero in the diagonal
            if (zero.isContained(augmentedMatrix[i][i][k])) {
                for (int j = i + 1; j < n; ++j) {
                    if (!zero.isContained(augmentedMatrix[j][i][k])) {
                        swapRows(augmentedMatrix, i, j);
                        break;
                    }
                }
            }

            // Scale the current row
            auto scale = augmentedMatrix[i][i][k];
            for (int j = 0; j < 2 * n; ++j) {
                augmentedMatrix[i][j][k] /= scale;
            }

            // Eliminate other rows
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    auto factor = augmentedMatrix[j][i][k];
                    for (int l = 0; l < 2 * n; ++l) {
                        augmentedMatrix[j][l][k] -= factor * augmentedMatrix[i][l][k];
                    }
                }
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    ComplexMatrix inverse;
    alloc(inverse, n, n, N);
    for (int k = 0; k < N; k++) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                inverse[i][j][k] = augmentedMatrix[i][j + n][k];
            }
        }
    }

    return inverse;
}